### functions to extract information form DP MCMC chains

##library to join arrays
library(abind)

#### auxiliary function to detect na values derived from possible errors ####
sum_nna<-function (a){sum( !is.na(a) ) }

#####  auxiliary function for numerical computation of integrals: Simpson method ####
simpson <- function(fun,h) {
  # numerical integral using Simpson's rule
  # assume a < b and n is an even positive integer
  n=length(fun)
  if (n == 2) {
    s <- fun[1] + 4*fun[2] +fun[3]
  } else {
    s <- fun[1] + fun[n] + 2*sum(fun[seq(3,n-2,by=2)]) + 4 *sum(fun[seq(2,n-1, by=2)] )
  }
  s <- s*h/3
  return(s)
}


#### Integrated absolute error function ####
IAE<-function(x,y,h){
  return( simpson(  abs(x-y),h ) )  
}

#### auxiliary function for numerical computation of integrals: trapezoidal integral ####
trapezoid_integral=function(v, h ){h*(mean(v[1],v[length(v)])+sum(v[2:(length(v)-1)]) )}

#### Bartlett–Priestley kernel ####
BPresker<-function(x, h  ){
  y=x  
  for(i in 1:length(x)){
    if(abs(x[i]) <h ){
      y[i]<-  (1- (x[i]/h)^2 )
    }else{
      y[i]<-   0
    }
  }
  return(y)  
}

##### Bartlett–Priestley kernel smoothing ####
kernBP<-function(x, Y, X, h){ 
  DD<-c()
  for(i in 1:length(x)){
    DD[i]<-    sum(   BPresker(x[i]-X   ,h)*Y  )          /sum(  BPresker(x[i]-X   ,h) ) 
  }
  return(DD)
}

#### Bartlett–Priestley kernel smoothing MSE ####
JBP=function(h, dataX, dataY){
  fhati<-rep(0,length(dataX))
  for(i in 1:length(dataX)){
    fhati[i]<- kernBP(x=dataX[i], X=dataX[-i], Y=dataY[-i], h = h)
  }
  return(mean( (fhati - dataY)^2 ))
}

#### function to extract the mode of components parameters ####
#this library is for clustering the parameter based on a Gaussian mixture
library(mclust) # loaded in the wrapper code BMARD_Replication.R





#auxiliar library to simulate AR(2) process
#library(TSA)

#### function to simualte AR2 mixture from the mixture components ####
# psi is a vector of phase parameter of the characteristic roots of the AR(2) 
# phi is a vector of weights of the mixture
# L is a vector of the log of the modulus of the characteristic roots of the AR(2) process
# n is the numer of points to simulate
# sd is the standard deviation of the time series
AR2mixsimulation<-function(psi, phi, L, n, sd ){
  #variables definition  
  phi1<- exp(-L)*2*cos(2*pi*psi)
  phi2<--exp(-2*L) 
  simindividual<-matrix(rep(1:n,length(psi)  ), nrow = length(psi) )
  
  for( i in 1:length(psi) ){
    simindividual[i,]<-sqrt(phi[i])*arima.sim(model=list(ar=c(phi1[i], phi2[i] ) )  , n=n)    
  }
  simmixture<-apply(simindividual, 2, sum)  
  
  return( list( individualmatrix=simindividual, simmixture=simmixture   )  )  
}


# The next functions extract the MCMC iterations that achieve high loglikelihood values,
# data is the list result from running the MBMARD function
# quant define the quantile level for the loglikelihood values .95 is used in our paper
# persamp define the percetage of the last MCMC iterations from the chains used in the computations before filtering by likelihood 
# qlow and qsup are the limits for the pointwise curves estimated for each chain and globally

## Computes the coherence PDC and the centers of the mixture components from the results of the MBMARD model
MBMARD_Coherence_PDC<-function( data , quant,percsamp,qlow, qsup ) {
  Nsamp= length( data[[1]]$Whittle )
  chains<-length(data)
  sampconsider<-round( Nsamp*percsamp,0)
  alllambda<-c()
  extweight<-c()
  psivecs<-c()
  BEWvecs<-c()  
  channelindicator<-c()
  N<-nrow(data[[1]]$Pgram)
  fourierdata<-data[[1]]$Pgram
  for(t in 1:chains ){
    if(t ==1){
      potentialhood<-c(data[[t]]$Whittle[(   sampconsider  +1 )  :Nsamp] )
    }else{
      potentialhood<-c(potentialhood, data[[t]]$Whittle[(   sampconsider  +1 )  :Nsamp] )  
    }
  }

  MLElist<- which( potentialhood > quantile(potentialhood,quant,na.rm = T  )  )
  
  select_by_chain<-c()
  index_by_chain<-c()
  myconst<-Nsamp-sampconsider
  for(t in 1:chains){
    
    alllambdatemp<-c()
    extweighttemp<-c()
    psivecstemp<-c()
    BEWvecstemp<-c()
  
    select_by_chain<- which( (myconst*(t-1)+1)<MLElist & MLElist<myconst*t  )
    if(length(select_by_chain)>0 ){
      index_by_chain<-    MLElist[select_by_chain] - myconst*(t-1)+ sampconsider
      curvesmatrixtemp<-array(data = NA,c(N,N, ncol(data[[1]]$Pgram ),   length(index_by_chain)) )
      GetC<-  GetCurvesBDP(Nsamp=length(index_by_chain), data[[1]]$omegas, psi_cur=data[[t]]$psi[index_by_chain,], BW_cur=data[[t]]$BW[index_by_chain,], Lambda_samps=data[[t]]$Lambda[index_by_chain], w=data[[t]]$weights[,,index_by_chain] )
        
      for(k in 1:length(index_by_chain) ){
        dftemp1<- data.frame( psi= data[[t]]$psi [index_by_chain[k],1:data[[t]]$Lambda[( index_by_chain[k])]], 
                              L= data[[t]]$BW[index_by_chain[k],1:data[[t]]$Lambda[( index_by_chain[k])] ] 
        )
        dftemp2<- cbind( dftemp1, t( data[[t]]$weights [, 1:data[[t]]$Lambda[( index_by_chain[k])  ] , index_by_chain[k] ]),data[[t]]$Lambda[( index_by_chain[k])  ]   )
        if(k==1){dftempall=dftemp2}else{
          dftempall=rbind(dftempall,dftemp2 )
        } 
        #compute curves
        curvesmatrixtemp[,,,k]<- GetC[[  k ]]  
        }
      
      
      alllambdatemp<-data[[t]]$Lambda[(index_by_chain)]
      alllambda<-c(alllambda, alllambdatemp)  
      extweight<-c(extweight,extweighttemp)
      psivecs<-c(psivecs,psivecstemp)
      BEWvecs<-c(BEWvecs,BEWvecstemp)
      #putting ell the data together
      if(!exists("curvesmatrix") ){curvesmatrix<-curvesmatrixtemp}else{
        curvesmatrix<-abind(curvesmatrix, curvesmatrixtemp )
      }    
    }# else{index_by_chain[[t]]=NA }
    
  } #end chains

  mediancurve=apply( curvesmatrix,c(1,2,3),median  )
  meancurve=apply( curvesmatrix,c(1,2,3),mean  )
  quantilcurve1<-apply( curvesmatrix,c(1,2,3),quantile,qsup  )
  quantilcurve2<-apply( curvesmatrix,c(1,2,3),quantile,qlow  )

  priormeans = as.numeric(names (   table(dftempall[,N+3] ))) 

  iclvec<-c()
  # compute the main components from potential number of clusters and maximize the integrated conditional likelihood
  for( g in  1: length( priormeans)  ){
    #      g=1

    k2 <- Mclust(dftempall[,1:(N+2)], G=priormeans[g]) 
    iclvec[g]<- k2$icl #to maximize
  }
  
  indexlambda <-which.max(iclvec)
  candidate_lambda<-priormeans[indexlambda]
  k2 <- Mclust(dftempall[,1:(N+2)],candidate_lambda) 
  centers<-k2$parameters$mean
  varcenters<-k2$parameters$variance$sigma
  mclustlhood<- k2$icl #to maximize
  sumweight<- apply(centers, 1,sum )
  for(h in 1:N){
    centers[(h+2), ]<-  centers[(h+2), ]/sumweight[(h+2)]   
  }  
   
  
  dfpdc<-dftempall[ which( dftempall[,(N+3)]==candidate_lambda ),]
  
limpdc <- dim(dfpdc)[1]/candidate_lambda

pdcarray<-array(data=0,  dim =c(N,candidate_lambda, length(data[[1]]$omegas) ,limpdc) )
coharray<-array(data=0,  dim =c(N,candidate_lambda, length(data[[1]]$omegas),limpdc ) )
for(l in 1:limpdc ){
  pdcarray[,,,l]<-PartialDirCoh(data[[1]]$omegas, t( dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,3:(N+2)] ), dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,1] , dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,2])  
  coharray[,,,l]<-spectralmatrix_all(data[[1]]$omegas, t( dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,3:(N+2)] ), dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,1] , dfpdc[((l-1)*candidate_lambda+1):(l*candidate_lambda)  ,2])$coherenceZX   
}  

medianpdc=apply( pdcarray,c(1,2,3),median  )
meanpdc=apply( pdcarray,c(1,2,3),mean  )
quantilpdc1<-apply( pdcarray,c(1,2,3),quantile,qsup  )
quantilpdc2<-apply( pdcarray,c(1,2,3),quantile,qlow  )

mediancoh=apply( coharray,c(1,2,3),median  )
meancoh=apply( coharray,c(1,2,3),mean  )
quantilcoh1<-apply( coharray,c(1,2,3),quantile,qsup  )
quantilcoh2<-apply( coharray,c(1,2,3),quantile,qlow  )

  return(list(
    mediancurve=mediancurve,
    meancurve=meancurve,
    lowcurve=quantilcurve1 ,
    supcurve=quantilcurve2,
    fourierdata=fourierdata,
    centers=centers,
    medianpdc=medianpdc,
    meanpdc=meanpdc,
    quantilpdc1=quantilpdc1,
    quantilpdc2=quantilpdc2,
    mediancoh=mediancoh,
    meancoh=meancoh,
    quantilcoh1=quantilcoh1,
    quantilcoh2=quantilcoh2
  ) )
}

#This function summarizes the centers and their variance of the mixture components from the results of the MBMARD model
MBMARD_components<-function( data , quant,percsamp,qlow, qsup ){
  
  Nsamp= length( data[[1]]$Whittle )
  chains<-length(data)
  sampconsider<-round( Nsamp*percsamp,0)
  alllambda<-c()
  extweight<-c()
  psivecs<-c()
  BEWvecs<-c()  
  channelindicator<-c()
  N<-nrow(data[[1]]$Pgram)
  fourierdata<-data[[1]]$Pgram
  for(t in 1:chains ){
    if(t ==1){
      potentialhood<-c(data[[t]]$Whittle[(   sampconsider  +1 )  :Nsamp] )
    }else{
      potentialhood<-c(potentialhood, data[[t]]$Whittle[(   sampconsider  +1 )  :Nsamp] )  
    }
  }
  
  MLElist<- which( potentialhood > quantile(potentialhood,quant,na.rm = T  )  )
  select_by_chain<-c()
  index_by_chain<-c()
  myconst<-Nsamp-sampconsider
  for(t in 1:chains){
    
    alllambdatemp<-c()
    extweighttemp<-c()
    psivecstemp<-c()
    BEWvecstemp<-c()
    select_by_chain<- which( (myconst*(t-1)+1)<MLElist & MLElist<myconst*t  )
    
    if(length(select_by_chain)>0 ){
      index_by_chain<-    MLElist[select_by_chain] - myconst*(t-1)+ sampconsider
      for(k in 1:length(index_by_chain) ){
        dftemp1<- data.frame( psi= data[[t]]$psi [index_by_chain[k],1:data[[t]]$Lambda[( index_by_chain[k])]], 
                              L= data[[t]]$BW[index_by_chain[k],1:data[[t]]$Lambda[( index_by_chain[k])] ] 
        )
        
        dftemp2<- cbind( dftemp1, t( data[[t]]$weights [, 1:data[[t]]$Lambda[( index_by_chain[k])  ] , index_by_chain[k] ]),data[[t]]$Lambda[( index_by_chain[k])  ]   )
        if(k==1){dftempall=dftemp2}else{
          dftempall=rbind(dftempall,dftemp2 )
        } 
        
      }
      
      alllambdatemp<-data[[t]]$Lambda[(index_by_chain)]
      alllambda<-c(alllambda, alllambdatemp)  
      extweight<-c(extweight,extweighttemp)
      psivecs<-c(psivecs,psivecstemp)
      BEWvecs<-c(BEWvecs,BEWvecstemp)
    }
    
  } #end chains
  
  priormeans = as.numeric(names (   table(dftempall[,N+3] )))
  
  toppmean<-max(priormeans)
  priormeans<- c(priormeans, c(1:10+toppmean ) )
  iclvec<-c()
  
  # compute the main components from potential number of clusters and maximize the integrated conditional likelihood
  for( g in  1: length( priormeans)  ){
    k2 <- Mclust(dftempall[,1:(N+2)], G=priormeans[g]) 
    iclvec[g]<- k2$icl #to maximize
  }
  
  indexlambda <-which.max(iclvec)
  candidate_lambda<-priormeans[indexlambda]
  k2 <- Mclust(dftempall[,1:(N+2)],candidate_lambda) 
  centers<-k2$parameters$mean
  varcenters<-k2$parameters$variance$sigma
  mclustlhood<- k2$icl #to maximize
  sumweight<- apply(centers, 1,sum )
  for(h in 1:N){
    centers[(h+2), ]<-  centers[(h+2), ]/sumweight[(h+2)]   
  }  
  
  return(list(
    centers=centers,
    varcenters=varcenters
  ) )
}

