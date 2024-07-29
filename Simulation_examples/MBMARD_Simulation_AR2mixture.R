# Multivariate Bayesian Mixture AutoRegressive Decomposition 
# EXAMPLE of a mixture of AR(2) processes
# by Guillermo C. Granados Garcia

#### Basic parameters ####

#number of mcmc iterations
Nsamp=1000 # note it is a small number just as a show case for a good estimation consider increase the mcmc samples
#time series length
size=1000
#sample frequency 
Sample_frequency=1000  # measured in Hertz (Hz), time points per second.

### loading and installing required packages

## next code to require the libraries is taken from the next article 
## https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("Rcpp", "TSA","stringr",
             "doParallel", "mclust","rstudioapi")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


####  sourcing the necessary code ####
codepath<-dirname(rstudioapi::getActiveDocumentContext()$path )
pathsplit= str_split_1(codepath,"/")
pathsplit= pathsplit[-length(pathsplit) ]
codepath=paste(pathsplit, collapse="/")
Rcpp::sourceCpp(paste(codepath, '/CPPcode/MuiltivariateBMARD.cpp',sep="" ) ) 
source(paste(codepath, '/AuxilliaryRcode/MBMARDchainssummaryCurvecomputation.R',sep="" ))



## The AR 2 mixture is the same as defined in Section 4

dimnodes<-7

psi<- (1/1000)*c(5, 30,60,  300)
L<-rep(.03, 4)
omegas =seq(1/size,.5-1/size ,by=1/size)

par1<- c(exp(-L[1])*2*cos(2*pi*psi[1]), -exp(-2*L[1])   )
par2<-c(exp(-L[2])*2*cos(2*pi*psi[2]), -exp(-2*L[2])   )
par3<-c(exp(-L[3])*2*cos(2*pi*psi[3]), -exp(-2*L[3])   )
par4<-c(exp(-L[4])*2*cos(2*pi*psi[4]), -exp(-2*L[4])   )
S1<-arima.sim(n = size, list(ar = par1 ), sd = 1, n.start = 10000)
S1<-(S1-mean(S1))/sd(S1) + rnorm(size, sd=.1)
S2<-arima.sim(n = size, list(ar = par2 ), sd = 1, n.start = 10000)
S2<-(S2-mean(S2))/sd(S2) + rnorm(size, sd=.1)
S3<-arima.sim(n = size, list(ar = par3 ), sd = 1, n.start = 10000)
S3<-(S3-mean(S3))/sd(S3)+ rnorm(size, sd=.1)
S4<-arima.sim(n = size, list(ar = par4 ), sd = 1, n.start = 10000)
S4<-(S4-mean(S4))/sd(S4)+ rnorm(size, sd=.1)
S5<-.4^.5*S1+.6^.5*S3
S5<-(S5-mean(S5))/sd(S5)+ rnorm(size, sd=.1)
S6<-.7^.5*S2+.3^.5*S4
S6<-(S6-mean(S6))/sd(S6)+ rnorm(size, sd=.1)
S7<-.3^.5*S3+.7^.5*S4
S7<-(S7-mean(S7))/sd(S7)+ rnorm(size, sd=.1)


allTS<-cbind(S1,S2,S3,S4,S5,S6,S7 )
FTallts <-t(mvfft(allTS)[1:499,]) /sqrt(size)

# optional:visualization of the spectral matrix raw data 
#matplot(abs(t(FTallts)), type="l")

#####  MBMARD parameters #####
#the frequencies at which the estiamtion is evaluated
omegas =seq(1/size,.5-1/size ,by=1/size)

# truncation of the Dirchlet process # of atoms
truncation=15

# initial number of components
m=10
# absolute squared value of the Fast Fourier matrix
periodmat<-abs(FTallts)^2
# initial location of the components we set in this example evenly across the interval (0,.5)
psiinit<-   seq(.01,.45,length.out=m)
#number of dimensions of the times series
n<-7
# initial Dirichlet process weights: in this example we set them randomly in (0,1)
Vinit<-array(runif(n*truncation ), dim=c(n,truncation )  )
# initial location of the Dirichlet process atoms: in this example we asign 1 per DP initial component
# aligned with the initial psi, the rest we asign them randomly
Zinit<-matrix(c(psiinit*2 , runif(truncation-m ) ), nrow =n, ncol =truncation, byrow = T )
# the metropolis hastings steps size for the DP atoms 
Linit<-runif(m,0.0001,.5)

# parameters of the MBMARD function 

#Pgram is a matrix containing the time series to be analyzed
#Nsamp is the number of posterior samples
#L is the truncation level for the BDP
#epsilon is a vector of uniform halfwidths for V and Z proposals
#epsilon_BW is a vector of uniform halfwidths for Bandwidth proposals length same as lambda_init
#SamplingRate is the sampling rate of the time series
#MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
#Sup is how often status updates are printed to the console, every Sup iterations
#Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
#Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
#either "expon", "poisson", or "flat" (only poisson is implemented so far)
#Lambda_max is the maximum number of Bernstein polynomial components to consider
#Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps (only one option is implemented so far)
#either "poisson" or "up1down1". If mispecifiec or "const" then lambda will be constant 
#alpha_const is a boolean to determine if the concetration parameter is set to constant
#alpha_init is the initial value for concentration parameter and the value for all samples if constant
#a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha is to be sampled
#q is the limit for the uniform distribution of the bandwidth parameters


MBMARD_test<- MBMARD( Tsize=size , omegas, Pgram=FTallts,periodmat,  Nsamp,  truncation,   epsilonV=.1,  epsilonZ=.05 ,  epsilon_BW=.1, Zinit=Zinit , Vinit=Vinit,
                         psi_init=psiinit,  BWinit=Linit,  
                         SamplingRate=1,  MaxFreq=0.5,  Sup=200, 
                         Lambda_init=m , Lambda_prior = "flat",
                         Lambda_max = 20, Lambda_prop_dist = "up1down1", alpha_const = F,
                         alpha_init = 1, a_alpha_L = .1, b_alpha_L = .1, q = 2.5,sigma=1, epsilon_sigma=.1
)




#monitor the loglikelihood
plot(MBMARD_test$Whittle, type="l")

#monitor the number of components
plot(MBMARD_test$Lambda)


##### extraction of main curves and components centers #####
## The Function accepts multiple chains as a list
## for this particular example only one chain is run so it must be insert into a list
## percsamp determines the burn-in period expressed as a percentage in (0,1)
MBMARD_test_list=list(MBMARD_test)
listmodes<- MBMARD_Coherence_PDC(data=MBMARD_test_list,quant=.9,percsamp=.5 ,qlow=.025, qsup=.975  )

# visualization of the spectral matrix estimation
#Mean curve = black curve
#lower credible interval curve = blue curve
#higher credible interval curve = green curve
#location of the components psi parameter estimated = violet vertical lines
#the location lines widths are determined by the weights of each component 
par(mfrow=c(n,n),mar=c(0,0,0,0) )
for(i in 1:n){
  for(j in 1:n){
    plot(omegas, listmodes$meancurve[i,j,],type = "l", lwd=2 , xaxt="n", yaxt="n")
    abline(v=listmodes$centers[1,], lwd=listmodes$centers[(j+2),]*4, col=6 )
  if(i==j){  
    lines(omegas, abs(  MBMARD_test$Pgram[i,]  )^2, col=2 , type="l" , lwd=1, lty=3 )  
}
    lines(omegas, listmodes$lowcurve[i,j,],type = "l", col=3, lwd=2  )
    lines(omegas, listmodes$supcurve[i,j,],type = "l",col=4, lwd=2  )
  
  }}


##### Partial Directed coherence curves #####
# here we plot the PDC from the latent components towards the observed series
# all plots are scaled in (0,1)

#number of all centers 
p=dim(listmodes$centers)[2]

par(mfrow=c(n,p),mar=c(0,0,0,0) )
for(i in 1:n){
  for(j in 1:p){
    plot(omegas, listmodes$meanpdc[i,j,],type = "l", lwd=1 , xaxt="n", yaxt="n",ylim = c(0,1) )
    lines(omegas, listmodes$quantilpdc1[i,j,],type = "l", col=2, lwd=2  )
    lines(omegas, listmodes$quantilpdc2[i,j,],type = "l",col=3, lwd=2  )
  }}





