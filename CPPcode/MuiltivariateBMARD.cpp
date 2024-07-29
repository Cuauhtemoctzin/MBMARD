#include <iostream>
#include <string>
#include <cmath>
#include <math.h> 
using namespace std;
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
#include <Rcpp.h>
using namespace Rcpp;

// Multivariate Bayesian Mixture AutoRegressive Decomposition (MBMARD)
// By Guillermo Cuauhtectzin Granados Garcia



//spectral spline function
// compute the mixture of AR(2) spectrums over a set of frequencies 
// [[Rcpp::export]]
vec IndividualSpline(const double& psi,const double& L, const vec& freq){
  //phi is a vector of component weights
  //psi is a vector of components AR2 phases
  //L is a vector of component log modulus each based on a AR2 process  
  double pi=3.14159265358979323846;
  vec fval(freq.size(),fill::zeros);
  vec Q(freq.size(),fill::zeros);
  for(int j =0; j< freq.size() ; j++ ){
    fval(j)+= 2.0*((((1-exp(-2.0*L))/(1+exp(-2.0*L)))*(  pow(1+exp(-2.0*L) ,2)  -4.0*(  pow( cos(2.0*pi*psi),2 )  )*exp(-2.0*L) ))/(  pow( 1-2.0*exp(-L)*cos(2.0*pi*psi)*cos(2.0*pi*freq[j])+exp(-2.0*L)*cos(4.0*pi*freq[j]), 2) +  pow (-2.0*exp(-L)*cos(2.0*pi*psi)*sin(2.0*pi*freq[j])+exp(-2.0*L)*sin(4.0*pi*freq[j]),2 )   ));  
  } // end for j
  return fval;
}

// compute the spectral matrix based on the splines and its parameters
// [[Rcpp::export]]
cube spectralmatrix(const vec& psi,const vec& L, const vec& freq){
  //psi is a vector of components AR2 phases
  //L is a vector of component log modulus each based on a AR2 process 
  // freq are the frequencies to compute the spectral matrix
  // assuming non-correlated components the matrix is diagonal in a 3d form
  int p= psi.size();
  cube specmat(p,p,freq.size(),fill::zeros ) ;
  vec temp(freq.size());
  for(int i=0 ;i< psi.size();i++ ){
    temp=IndividualSpline( psi(i), L(i),  freq);
    specmat.tube(i,i)=temp; 
  }
  return(specmat);
}

// central cube of the nodes after mixing the latent sources through the weight matrix 
// [[Rcpp::export]]
cube spectralmatrix_nodes( const vec& freq, mat weights, cube specmat){
  //   m is the number of nodes
  int m=weights.n_rows;
  mat f1(m,m);
  mat f2(m,specmat.n_cols);
  vec a(specmat.n_cols);
  cube specnodes(m,m,freq.size(),fill::zeros )  ;  
  for(int i =0;i<freq.size();i++){
  a=specmat.slice(i).diag();
    for(int j=0;j<m;j++ ){
    f2.row(j)= weights.row(j)%a.t();
    }
    f1=f2* weights.t();
    specnodes.slice(i)=  f1;
  }
  return specnodes;
}


// csepctral matrix considering all the sources, the observed data and the latent components 
// [[Rcpp::export]]
List spectralmatrix_all( const vec& freq, mat weights, const vec& psi,const vec& L){
  //   m is the number of nodes
  // the coherence is the lower triangular part of the matrix of coherence 
  // but set all in a set of columns 
  // this structure match other coherence functions in R
  // the structure for the conherence between Z and X is a cube 
  
  cube specmat=  spectralmatrix( psi, L, freq);
  int m=weights.n_rows;
  int n=weights.n_cols;
  mat f2(m,specmat.n_cols);
  mat f1(m,m);
  vec a(specmat.n_cols);
  cube spec_all(m+n,n+m,freq.size(),fill::zeros )  ;
  mat A(m,n);
  int dim1=m*(m-1)/2;
  mat cohmatXX(freq.size(),dim1 ) ;
  cube cohmatZX(m,n,freq.size() ) ;
  for(int i =0;i<freq.size();i++){
    a=specmat.slice(i).diag();
    for(int j=0;j<m;j++ ){
      f2.row(j)= weights.row(j)%a.t();
    }
    f1=f2* weights.t(); 
    spec_all.subcube(0,0,i,n-1,n-1,i)= specmat.slice(i)  ;  
    spec_all.subcube(n,n,i,m+n-1,m+n-1,i)= f1; 
    A=weights*specmat.slice(i);
    
    int flag=0;
    for(int k =1;k<m;k++){
      for(int j =0;j<k;j++){
        cohmatXX(i,flag)= f1(k,j)/ sqrt(f1(k,k)*f1(j,j))  ;
        flag++;
      }
    }
    for(int k =0;k<m;k++){
      for(int j =0;j<n;j++){
        cohmatZX(k,j,i)= A(k,j)/ sqrt(specmat(j,j,i ) * f1(k,k))  ;
      }
    }
  // next is z x  

    spec_all.subcube(n,0,i,m+n-1,n-1,i)= A  ;
    spec_all.subcube(0,n,i,n-1,m+n-1,i)= A.t()  ;
  }
  //return spec_all;
  return List::create(
    _["spec_all"] = spec_all,// _["inveigen"] = inveigen);
    _["coherenceXX"] = cohmatXX,
    _["coherenceZX"] = cohmatZX) ;
  
  }


//function to compute the module of the coefficients
// [[Rcpp::export]]
mat ModulePardirVec( const vec& freq, vec psi, vec bw){
  //\Phi_j(\omega)=\phi_{j1}e^{-i2\pi\omega}+ \phi_{j2}e^{-i4\pi\omega}  
  //|\pi_{ij}(\omega)|= \frac{\lambda_{ij} |\Phi_j(\omega)| }{\sqrt{|\Phi_j(\omega) |^2 \sum_{i=1}^n \lambda^2_{ij} +|1-\Phi_j(\omega)|^2  } }  
  double pi=3.14159265358979323846;
int  kf= freq.size();
int comp=psi.size();
mat phimodule(kf,comp);   
vec phi1(comp);  
vec phi2(comp);  
for(int i=0;i<comp;i++){
  phi1(i)= exp(-bw(i))*2*cos(2*pi*psi(i)) ;
  phi2(i)= -exp(-2*bw(i)) ;
 for(int j=0;j<kf;j++){
   phimodule(j,i)=  pow(phi1(i)*cos(2*pi*freq(j))+ phi2(i)*cos(4*pi*freq(j)),2)  + pow(phi1(i)*sin(2*pi*freq(j))+ phi2(i)*sin(4*pi*freq(j)),2);
 }
}
return(phimodule);
}

//function to copmute the other module of the pdc
// [[Rcpp::export]]
mat Moduleoneminus( const vec& freq, vec psi, vec bw){
  //\Phi_j(\omega)=\phi_{j1}e^{-i2\pi\omega}+ \phi_{j2}e^{-i4\pi\omega}  
  //|\pi_{ij}(\omega)|= \frac{\lambda_{ij} |\Phi_j(\omega)| }{\sqrt{|\Phi_j(\omega) |^2 \sum_{i=1}^n \lambda^2_{ij} +|1-\Phi_j(\omega)|^2  } }  
  double pi=3.14159265358979323846;
  int  kf= freq.size();
  int comp=psi.size();
  mat phimodule(kf,comp);   
  vec phi1(comp);  
  vec phi2(comp);  
  for(int i=0;i<comp;i++){
    phi1(i)= exp(-bw(i))*2*cos(2*pi*psi(i)) ;
    phi2(i)= -exp(-2*bw(i)) ;
    for(int j=0;j<kf;j++){
      phimodule(j,i)=  pow(1-phi1(i)*cos(2*pi*freq(j))- phi2(i)*cos(4*pi*freq(j)),2)
      + pow(phi1(i)*sin(2*pi*freq(j))+ phi2(i)*sin(4*pi*freq(j)),2);
    }
  }
  return(phimodule);
}


//function to compute the partial directed coherence of the MBMARD model 
// [[Rcpp::export]]
cube PartialDirCoh( const vec& freq, mat weights,  vec psi, vec bw){
int N = weights.n_rows;
int comp = weights.n_cols;
int nfq = freq.size();  
cube PDC(N,comp, nfq);
vec sumsq(comp);
///vector sum of squares of the weight matrix
mat phimodule=ModulePardirVec( freq,  psi,  bw) ;
mat phimodonemin=Moduleoneminus( freq,  psi,  bw) ;
mat sqweight= weights%weights;
//b has the same size of the components
vec b = sum(sqweight,0).t();
for(int i =0;i<N;i++  ){
  for(int j =0;j<comp;j++  ){
    for(int k =0;k<nfq;k++  ){
     PDC(i,j,k)= weights(i,j)*sqrt(phimodule(k,j) )/sqrt( phimodule(k,j)*b(j) +phimodonemin(k,j)   ) ;
    }
  }
}
  return(PDC);
}



//matrix breaking representation 
//givven beta V build wieghts PI_mh m nodes , h trunctation index 
// [[Rcpp::export]]
mat GetWeightsMult( mat Vmat){
  //V is the same dimension  m nodes  h truncation(atoms)
  int L=Vmat.n_cols ;
  vec np(1);
  np(0)=0;
  mat Pmat(Vmat.n_rows,Vmat.n_cols,fill::zeros);
  rowvec V1(L);
  vec V(L);
  vec logV(L);
  vec aux(L-1);
  vec CumProds(L);
  vec out(L);
  for(int i =0;i<Vmat.n_rows;i++){
    V1=Vmat.row(i);
    V=V1.t();
    logV=log(V);
    aux=  V.subvec(0,L-2) ;
    CumProds= join_cols(np, cumsum( log(1-aux) ) );
    out= exp(CumProds+logV) ;
    Pmat.row(i)= out.t();
  }
  return Pmat;
}

//Function to calculate BDP weights from stick break proportions V for truncation level L for a matrix results a matrix NXP
// [[Rcpp::export]]
mat GetPhi(const mat& Z,const mat& Pmat, int Lambda_now){
  //Z is matrix L x M of atoms and the number of atoms determine the weights to sum 
  //P is a matrix with individual weights for atoms for factors dimension p and nodes dimensions m
  //L is the level of truncation of the DP 
  //Lambda_now is the current number of latent factors
  int L=Pmat.n_cols;
  mat Phimat(Pmat.n_rows,Lambda_now,fill::zeros );
  rowvec P1;
  vec P;
  rowvec Z1;
  vec Z2;
  vec Phi(Lambda_now);
  for(int i=0;i<Pmat.n_rows;i++ ){
    P1 = Pmat.row(i);
    Z1= Z.row(i) ;
    Z2= Z1.t() ;
    P = P1.t();
    Phi.fill(0);
    for(int l=0; l<L; l++){
      Phi[ floor(Z2[l]*Lambda_now)  ] += P[l];
    }
    Phimat.row(i)=Phi.t();
  }//end i
  return Phimat;
}


//Function to calculate BDP weights from stick break proportions V for truncation level L for a matrix results a matrix NXP
// this function in particular aligns the weights with the components of the partition 
// [[Rcpp::export]]
mat GetPhi_part(const mat& Z,const mat& Pmat, int Lambda_now, vec partition){
  //Z is matrix L x M of atoms and the number of atoms determine the weights to sum 
  //P is a matrix with individual weights for atoms for factors dimension p and nodes dimensions m
  //L is the level of truncation of the DP 
  //Lambda_now is the current number of latent factors
  int L=Pmat.n_cols;
  mat Phimat(Pmat.n_rows,Lambda_now,fill::zeros );
  rowvec P1;
  vec P;
  rowvec Z1;
  double phitemp;
  vec Z2;
  vec Phi(Lambda_now);
  for(int i=0;i<Pmat.n_rows;i++ ){
    P1 = Pmat.row(i);
    Z1= Z.row(i) ;
    Z2= Z1.t() ;
    P = P1.t();
    Phi.fill(0);
    for(int k =0; k<Lambda_now;k++ ){  
      phitemp=0;
      for(int l=0; l<L; l++){
        if( partition[k]*2 <= Z2[l] &  Z2[l]<= partition[k+1]*2 ){
          phitemp += P[l];}
      }
      Phi[k]=  phitemp;
    }
    Phimat.row(i)=Phi.t();
  }//end i
  return Phimat;
}





//Function to propose new location parameter (phase of AR(2)) for a subinterval j with k total components 
// [[Rcpp::export]]
vec Newpsij(const vec& psi ,   const double& epsilonpsi, int j, vec partition ){
  //current phi vector 
  //current psi vector 
  //current L vector
  // j is the index for the subinterval to sample the location parameter psi[j]
  // k total of components considered
  //epsilon is a vector of uniform half widths for proposals chose considering the size of the subinterval that is 1/k
  
  double prop=psi(j);
  vec propvec=psi ;
  double psinew;
  double v= (partition(j+1)-partition(j))/2.01    ; //1.0/static_cast<double>(2*k); 
  //double d=static_cast<double>(j)/static_cast<double>(2*k);
  //double z= static_cast<double>(j-1)/static_cast<double>(2*k);
  //propose new value
  prop = prop +v*(2*randu()-1);
  //take modulus inside the subinterval to ensure it falls inside the subinterval specified
  if(prop < partition(j) ){psinew=v+prop;}else{if(prop < partition(j+1)){psinew=prop;}else{psinew=prop-v;}}
  propvec(j)=psinew ;
  //propvec[j]= ;
  return propvec;
}


//Function to propose new bandwidth (log-modulus of AR(2) roots) for a subinterval j with k total components 
// [[Rcpp::export]]
vec NewL(const vec& L,const double& epsilon,double dim, double q){
  //VZ is vector of current weights or atoms
  //epsilon is a numerical value of uniform half widths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  
  vec Lnew=L;
  //propose new value
  double prop=L(dim)+epsilon*(2*randu()-1);
  //take modulus 1 to ensure it falls between 0 and 1
  if(prop < 0){Lnew(dim)=q+prop;}else{if(prop < q){Lnew(dim)=prop;}else{Lnew(dim)=prop-q;}}
  return Lnew;
}




//whittle log likelihood for multivariate time series using the Sherman-morrison-woodbury approach to reduce computational costs

// [[Rcpp::export]]
double WhittleSHMOWOOD(cx_mat fxx, mat weights, const vec& psi,const vec& L, const vec& freq, double sigma  ){
 
 double logdet=0;
  // sigma denotes the variance of the noise component
  cube specmat=spectralmatrix(psi, L,  freq);  
  int nchannels= fxx.n_rows;  
// Sherman-morrison-woodbury approach 
// (A+XBXT)^-1 =A^-1 - A^-1 X(B^-1 +XTA^-1X )^-1XTA^1

// A is the matrix with parameter sigma its inverse is diagonal with values 1/sigma
vec diagsigma(nchannels);  
diagsigma.fill(1/sigma) ;
mat Siginv=  diagmat(diagsigma);
// B is the matrix with the diagonal spectral matrix of AR(2) then the inverse is a diagonal cube
cube specmatinv(specmat.n_rows, specmat.n_cols, specmat.n_slices ); 
vec diagAR2(specmat.n_rows); 
vec invdiagAR2(specmat.n_rows);
for(int i =0;i<freq.size();i++){
  diagAR2=specmat.slice(i).diag();
  for(int j =0;j<specmat.n_rows;j++){
    invdiagAR2(j)=1/diagAR2(j) ;
  }
  specmatinv.slice(i)=diagmat(invdiagAR2) ;
}

// is important to compute the inverse and the determinant of this matrix X^T A^-1X
mat impmat= weights.t()*Siginv* weights  ;
mat sumimpmat(specmatinv.n_rows,specmatinv.n_cols);
  double fk=0;
  mat f2(Siginv.n_rows,Siginv.n_cols);
  cx_vec v1(fxx.n_rows);
  cx_vec v2(fxx.n_rows);
  for(int i =0;i<freq.size();i++){
// new procedure by using the sh-mor-wood formula
sumimpmat=specmatinv.slice(i) + impmat ;
f2=  Siginv- Siginv*weights * inv( sumimpmat  ) * weights.t() * Siginv;
 logdet =logdet+ sum(  log( specmat.slice(i).diag() ) ) +  log( det( sumimpmat  ) ) ;     
    v1= fxx.col(i);
    v2= v1.t()*f2*v1;
    fk=fk-  v2(0).real() ;
    //  }
  }
  fk=fk-logdet-  freq.size()*fxx.n_rows*log( sigma )   ;
  return fk;
  // note that is necessary to add the priors over all the parameters specially K the number of latent processes
  //  check if its possible to add the integrated likelihood, as a penalization for the mixture model  
}


//this function penalized the sum of squares of the diagonals (periodograms and the estimated)
// alow to have better fit towards the diagonals
// [[Rcpp::export]]
double Whittle_single_penalized(cx_mat fxx,mat periodmat , mat Vmat,const mat& Z, vec partition, const vec& psi,const vec& L, const vec& freq, double sigma  ){
  // the weight matrix is computed from the DP weights 
  /// periodmat has the smae dimension as fxx but are the individual periodogram values
  mat Pmat=GetWeightsMult( Vmat);
  int Lambda_now=psi.size() ;
  mat weights=  GetPhi_part( Z, Pmat,  Lambda_now,  partition);
  double logdet=0;
  // sigma denotes the variance of the noise component
  cube specmat=spectralmatrix(psi, L,  freq);  
  int nchannels= fxx.n_rows;  
  vec diagsigma(nchannels);  
  diagsigma.fill(1/sigma) ;
  mat Siginv=  diagmat(diagsigma);
  //this matrix XTA^-1X   si importantto compute the inverse and the determinant
  cube specmatinv(specmat.n_rows, specmat.n_cols, specmat.n_slices ); 
  vec diagAR2(specmat.n_rows); 
  vec invdiagAR2(specmat.n_rows);
    mat impmat= weights.t()*Siginv* weights  ;
  mat sumimpmat(specmatinv.n_rows,specmatinv.n_cols);
  cube specmatnodes=spectralmatrix_nodes(freq,weights,specmat);
  double penalization=0;
  double fk=0;
  mat f2(Siginv.n_rows,Siginv.n_cols);
  cx_vec v1(fxx.n_rows);
  cx_vec v2(fxx.n_rows);
  // B is the matrix with the diagonal spectral matrix of AR(2) then the inverse is a diagonal cube
  for(int i =0;i<freq.size();i++){
    diagAR2=specmat.slice(i).diag();
    for(int j =0;j<specmat.n_rows;j++){
      invdiagAR2(j)=1/diagAR2(j) ;
    }
    specmatinv.slice(i)=diagmat(invdiagAR2) ;
  }
  for(int i =0;i<freq.size();i++){
    // new procedure by using the sh-mor-wood formula
    sumimpmat=specmatinv.slice(i) + impmat ;
    f2=  Siginv- Siginv*weights * inv( sumimpmat  ) * weights.t() * Siginv;
    logdet =logdet+ sum(  log( specmat.slice(i).diag() ) ) +  log( det( sumimpmat  ) ) ;     
    v1= fxx.col(i);
    v2= v1.t()*f2*v1;
    fk=fk-  v2(0).real() ;
  }
  for(int j =0;j<nchannels;j++){
    vec tempvec=specmatnodes.tube(j,j);
    penalization=penalization +   sum( pow(periodmat.row(j).t()  -tempvec  ,2 ) ) ;
  }
  fk=fk-logdet-  freq.size()*fxx.n_rows*log( sigma )  -penalization ;
  return fk;
}


//function to create a new partition points depending if the size grow by 1 or decrease by 1 (many cases!!)
// [[Rcpp::export]]
List newpartition(vec partition, int growth, int i, vec BW, vec psi, mat V,const mat& Z,const double& epsBW, double q ){
  vec newpart ;
  vec psiinit;
  vec BWinit;
  mat P=GetWeightsMult(V);
  double explorepsi=2.1;
  int N =partition.n_elem;
  if(N==2){i=0;}
  if(i==N-2&&growth!=1&&N!=2){i=N-3;}
  int Lambda_now=N-1;
  double b= partition(i+1);
  double a= partition(i);
  mat phimat= GetPhi_part( Z, P, Lambda_now, partition);
  
  vec phi(Lambda_now);
  for(int i=0;i<Lambda_now;i++ ){
    phi(i)= sum(phimat.col(i) )/static_cast<double>(phimat.n_rows) ; 
  }
  
  double newpoint=  randu()*(b-a )+a   ;
  vec np(1);
  vec inter(1);
  np(0)=newpoint;
  if (growth==1){//birth
    newpart=join_cols(partition.subvec(0,i), np);
    newpart=join_cols(newpart,  partition.subvec(i+1,N-1));
    if(N!=2){
      if(psi(i)<newpoint ){
        if(i==N-2){
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols( psi, inter);
          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
          inter(0)=BW(i) ;
          BWinit=join_cols(BW ,inter);
          BW=NewL(BWinit,   epsBW  ,i ,q);
        }else{
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i+1,N-2));
          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
          inter(0)=BW(i) ;
          BWinit=join_cols(BW.subvec(0,i),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i+1,N-2));
          BW=NewL(BWinit,   epsBW  ,i, q);
        }  
      }else{// case psi(i)>=newpoint 
        if(i==0){
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   inter;
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i))/explorepsi , i, newpart );
          inter(0)=BW(i) ;
          BWinit= inter;
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2)  );
          BW=NewL(BWinit,   epsBW  ,i , q);
        }else{
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i-1),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i) )/explorepsi , i, newpart );
          inter(0)=BW(i) ;
          BWinit=join_cols(BW.subvec(0,i-1),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2));
          BW=NewL(BWinit,   epsBW  ,i , q);
        }
      }
    }else{// case N=2
      if(psi(i)<newpoint ){
        
        inter(0)=(newpoint +partition(i+1))/2.0 ;
        psiinit=   join_cols(psi,inter);
        psi=Newpsij( psiinit ,  (newpart(i+2)-newpart(i+1))/explorepsi , i+1, newpart );
        inter(0)=BW(i) ;
        BWinit=join_cols(BW,inter);
        BW=NewL(BWinit,   epsBW  ,i+1 , q);        
      }else{
        inter(0)=(newpoint +partition(i))/2.0 ;
        psiinit=   join_cols(inter,psi);
        psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
        inter(0)=BW(i) ;
        BWinit=join_cols(inter, BW);
        BW=NewL(BWinit,   epsBW  ,i , q);
        
      }
    }
  }else{if(N!=2){//death
    if (i==0){
      np(0)=partition(0);
      newpart =join_cols<mat>( np,partition.subvec( 2,N-1) );
      
      
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)))*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1)   ;
      }
      if(N==3){
        //   psi = inter;
        psi=  Newpsij( inter ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
        
      }else{
        
        psi =join_cols<mat>( inter,psi.subvec( 2,N-2) );
        psi=  Newpsij( psi ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
      }
      inter(0)= min(BW(i),BW(i+1) )   ;
      if(N==3){
        BW=inter;
        BW=NewL(inter, epsBW, i , q);
        
        
      }else{
        BW =join_cols<mat>( inter,BW.subvec( 2,N-2) );
        BW=NewL(BW, epsBW, i , q);
      }
    }else{if(i==N-3){
      np(0)=partition(N-1);
      newpart=join_cols<mat>(  partition.subvec(0,N-3), np  );
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)) )*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1);}
      
      psi =join_cols<mat>( psi.subvec( 0,N-4),inter  );
      psi=Newpsij(psi ,  (newpart(i+1)-newpart(i))/explorepsi , i,  newpart);
      inter(0)= min(BW(i),BW(i+1) );
      
      
      BW =join_cols<mat>( BW.subvec( 0,N-4), inter );
      BW=NewL(BW, epsBW, i , q);
      
    }else{
      newpart =join_cols<mat>( partition.subvec(0,i),partition.subvec( i+2,N-1) );
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)))*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1)   ;}
      psiinit =join_cols<mat>( psi.subvec( 0,i-1),inter  );
      psi =join_cols<mat>(psiinit, psi.subvec( i+2,N-2)  );
      psi=Newpsij(psi ,  (newpart(i+1)-newpart(i))/explorepsi , i,  newpart);
      
      inter(0)= min(BW(i),BW(i+1) )   ;
      BWinit =join_cols<mat>(BW.subvec( 0,i-1), inter );
      BW =join_cols<mat>(BWinit, BW.subvec( i+2,N-2) );
      BW=NewL(BW, epsBW, i , q);
      
    }
    }
  }else{newpart=partition;}
  }
  return List::create(
    _["partition"] = newpart,
    _["psi"] = psi,
    _["BW"] = BW);
}





//function to create a new partition points depending if the size grow by 1 or decrease by 1 
// [[Rcpp::export]]
List newpartitionOPT(vec partition, int growth, int i, vec BW, vec psi, mat V,const mat& Z,const double& epsBW, double q, cx_mat fxx ,  mat periodmat ,  vec freq, double sigma ){
  double factorBW=.5;
  vec newpart ;
  vec psiinit;
  vec BWinit;
  mat P=GetWeightsMult(V);
  double explorepsi=2.1;
  int N =partition.n_elem;
  if(N==2){i=0;}
  if(i==N-2&&growth!=1&&N!=2){i=N-3;}
  int Lambda_now=N-1;
  double b= partition(i+1);
  double a= partition(i);
  mat phimat= GetPhi_part( Z, P, Lambda_now, partition);
  
  vec phi(Lambda_now);
  for(int i=0;i<Lambda_now;i++ ){
    phi(i)= sum(phimat.col(i) )/static_cast<double>(phimat.n_rows) ; 
  }
  
  double newpoint=  randu()*(b-a )+a   ;
  vec np(1);
  vec inter(1);
  np(0)=newpoint;
  if (growth==1){//birth
    newpart=join_cols(partition.subvec(0,i), np);
    newpart=join_cols(newpart,  partition.subvec(i+1,N-1));
    if(N!=2){
      if(psi(i)<newpoint ){
        if(i==N-2){
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols( psi, inter);

          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
          inter(0)=randu()*factorBW ;
          BWinit=join_cols(BW ,inter);
          BW=NewL(BWinit,   epsBW  ,i ,q);
        }else{
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i+1,N-2));
          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
          inter(0)=randu()*factorBW ;
          BWinit=join_cols(BW.subvec(0,i),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i+1,N-2));
          BW=NewL(BWinit,   epsBW  ,i, q);
        }  

      }else{// case psi(i)>=newpoint 
        if(i==0){
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   inter;
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i))/explorepsi , i, newpart );
          inter(0)=randu()*factorBW ;
          BWinit= inter;
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2)  );
          BW=NewL(BWinit,   epsBW  ,i , q);
        }else{
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i-1),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i) )/explorepsi , i, newpart );
          inter(0)=randu()*factorBW ;
          BWinit=join_cols(BW.subvec(0,i-1),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2));
          BW=NewL(BWinit,   epsBW  ,i , q);
        }
      }
    }else{// case N=2
      if(psi(i)<newpoint ){
        
        inter(0)=(newpoint +partition(i+1))/2.0 ;
        psiinit=   join_cols(psi,inter);
        psi=Newpsij( psiinit ,  (newpart(i+2)-newpart(i+1))/explorepsi , i+1, newpart );
        inter(0)=randu()*factorBW ;
        BWinit=join_cols(BW,inter);
        BW=NewL(BWinit,   epsBW  ,i+1 , q);        
      }else{
        inter(0)=(newpoint +partition(i))/2.0 ;
        psiinit=   join_cols(inter,psi);
        psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
        inter(0)=randu()*factorBW ;
        BWinit=join_cols(inter, BW);
        BW=NewL(BWinit,   epsBW  ,i , q);
        
      }
      
      
      
    }
  }else{if(N!=2){//death
    if (i==0){
      np(0)=partition(0);
      newpart =join_cols<mat>( np,partition.subvec( 2,N-1) );
      // choose one of the two psi
      
      //  compute likelihood 
  
      //  first case
      inter=psi(0) ;
      vec psi1 =join_cols<mat>( inter,psi.subvec( 2,N-2) );
      inter=BW(0) ;
      vec BW1 =join_cols<mat>( inter,BW.subvec( 2,N-2) );
      
      //likelihood
      mat phimat1= GetPhi_part( Z, P, Lambda_now-1, newpart);
      //newpart
      double lhood1 = Whittle_single_penalized( fxx, periodmat,V, Z, newpart , psi1 , BW1, freq, sigma);
      
      vec psi2 =psi.subvec( 1,N-2) ;
      vec BW2 = BW.subvec( 1,N-2) ;
      double lhood2 = Whittle_single_penalized( fxx, periodmat,V, Z, newpart   , psi2 , BW2, freq, sigma);
      
      if(lhood1<lhood2){ psi=psi2;BW=BW2; }else{
        psi=psi1; BW=BW1;
      }
      
      
      
      //    if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
      //      inter(0)= (phi(i)/(phi(i) + phi(i+1)))*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1)   ;
      //    }
      if(N==3){
        //   psi = inter;
        psi=  Newpsij( inter ,  (newpart(i+1)-newpart(i))/explorepsi , i, newpart );
        
      }
      
    }else{if(i==N-3){
      np(0)=partition(N-1);
      newpart=join_cols<mat>(  partition.subvec(0,N-3), np  );
      
      inter=psi(N-2) ;
      vec psi1 =join_cols<mat>( psi.subvec( 0,N-4), inter );
      inter=BW(N-2) ;
      vec BW1 =join_cols<mat>( BW.subvec( 0,N-4),inter );
      
      //likelihood
      mat phimat1= GetPhi_part( Z, P, Lambda_now-1, newpart);
      double lhood1 = Whittle_single_penalized( fxx,  periodmat,V, Z, newpart , psi1 , BW1, freq, sigma);
      
      vec psi2 =psi.subvec( 0,N-3) ;
      vec BW2 = BW.subvec( 0,N-3) ;
      
      double lhood2 = Whittle_single_penalized( fxx,  periodmat,V, Z, newpart , psi2 , BW2, freq, sigma);
      
      if(lhood1<lhood2){ psi=psi2;BW=BW2; }else{
        psi=psi1; BW=BW1;
      }    
      
      
      
    }else{
      newpart =join_cols<mat>( partition.subvec(0,i),partition.subvec( i+2,N-1) );
      
      vec psi1 =join_cols<mat>( psi.subvec( 0, i  ), psi.subvec( i+2,N-2) );
      
      vec BW1 =join_cols<mat>( BW.subvec( 0,i),BW.subvec( i+2,N-2) );
      
      //likelihood
      mat phimat1= GetPhi_part( Z, P, Lambda_now-1, newpart);
      double lhood1 = Whittle_single_penalized( fxx,  periodmat,V, Z, newpart , psi1 , BW1, freq, sigma);
      
      vec psi2 =join_cols<mat>( psi.subvec( 0, i-1  ), psi.subvec( i+1,N-2) ) ;
      vec BW2 = join_cols<mat>( BW.subvec( 0,i-1),BW.subvec( i+1,N-2) );
      
      double lhood2 = Whittle_single_penalized( fxx,  periodmat,V, Z, newpart , psi2 , BW2, freq, sigma);
      
      if(lhood1<lhood2){ psi=psi2;BW=BW2; }else{
        psi=psi1; BW=BW1;
      }
      
    }
    }
  }else{newpart=partition;}
  }
  return List::create(
    _["partition"] = newpart,
    _["psi"] = psi,
    _["BW"] = BW);
}







//function to get the width  of a subinterval 
//used to compute uniform densities on the M-H steps
// [[Rcpp::export]]
double Getwidth(vec partition , int i){
  double width= partition[i+1]- partition[i];
  return  width;
}

// replicate the function "sample" from R
// [[Rcpp::export]]
NumericVector mysample(NumericVector integervec, NumericVector probs){
  // calling rnorm()
  Function f("sample");   
  
  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(Named("x")= integervec, Named("size")=1, Named("prob")=probs);
}

//random selection of the component in death tend to join small intervals and in growth tend to split big intervals
//in a death-birth process
// [[Rcpp::export]]
vec rand_selection(vec partition, int K){
  //K is the growth death flag: is 1 when to growth 1 subinterval
  //creation of integer
  NumericVector widthvec(partition.size()-1);
  NumericVector integervec(partition.size()-1);
  for(int i=0;i<partition.size()-1;i++){
    integervec[i]= i ;
  }
  vec widthvectemp(partition.size()-1);
  NumericVector ret(1);
  double mysum;
  if(K==1){  // growth: 
    for(int i=0;i<partition.size()-1;i++){
      widthvectemp[i]= partition[i+1]-partition[i] ;
    }
    mysum=sum(widthvectemp) ;
    for(int i=0;i<partition.size()-1;i++){
      widthvec[i]=widthvectemp[i]/mysum;
    }
    ret =  mysample(integervec, widthvec );
  } else{
    for(int i=0;i<partition.size()-1;i++){
      widthvectemp(i)= 1.0/(partition(i+1)-partition(i)) ;
    }
    mysum=sum(widthvectemp) ;
    for(int i=0;i<partition.size()-1;i++){
      widthvec(i)=widthvectemp(i)/mysum;
    }
    ret = mysample(integervec ,widthvec);
  } 
  return ret;  
}



//function to sample number of components, Lambda, using MH step
// in a death-birth process
// [[Rcpp::export]]
List sample_Lambda(int N, int L, const vec& omegas,const vec& psi_cur,const vec&  BW_cur, const cx_mat& Pgrams, mat periodmat, const mat& P_cur, const mat& V_cur,
                   const mat& Z_cur, const double& epsBW, int Lambda_cur, string Lambda_prior,  
                   int Lambda_max, double MaxFreq,  vec partition, double q, double sigma, double lhood ){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each frequency
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column for each subject
  //Lambda_cur is the current number of Bernstein polynomial components
  //Lambda_prior is the prior for the number of Bernstein polynomial components
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //MaxFreq is the maximum frequency to consider for the spectral density
  //Lambda_prop_dist is the proposal distribution to be used for Lambda, up1down1 or poisson. If incorrectly/not specified Lambda will be constant
  //DBDP is a boolean to determine if DP atoms are dependent n covariates
  // i is the subinterval chosen to make the partition or the death 
 // int M=partition.n_elem-2;
 //likehood is the value of the whittleloglikelihood
 
  List  newpart;
  int M=partition.n_elem-2;
  int Lambda_new;
  int Lambda_prop = Lambda_cur;
  double accept_Lambda = 0;
  int growth;
  int i;
  vec   newp=partition;
  vec   newpsi=psi_cur;
  vec   newBW=BW_cur;

  if(Lambda_cur==2){Lambda_prop=3;
    growth=1;
    i=0;
    newpart= newpartitionOPT( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW  ,q,  Pgrams,periodmat,  omegas,  sigma);
  //  newpart= newpartition( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW  ,q);
    
    vec psinew=newpart[1];
    vec np=newpart[0];
  }
  else{
    if(Lambda_cur==Lambda_max){Lambda_prop=Lambda_max-1;
      growth=0;
      // mandatory death
     //  i=rand_selection(partition, growth)(0);
      i=randi<int>( distr_param(0, M));
      newpart= newpartitionOPT( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW  ,q,  Pgrams,periodmat,  omegas,  sigma);
      vec psinew=newpart[1];
      vec np=newpart[0];
    }else{
      if(randu()<0.5){Lambda_prop=Lambda_cur-1;
        growth=0;
       
        ///DEATH change to a better selection INVERSE proportional to parittion widths  join small subinterval
        //  i=rand_selection(partition, growth)(0);
        i=randi<int>( distr_param(0, M));
        newpart= newpartitionOPT( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW  ,q,  Pgrams,periodmat,  omegas,  sigma);
        
        vec psinew=newpart[1];
        vec np=newpart[0];
        
      }else{Lambda_prop=Lambda_cur+1;
        growth=1;
        
        ///GROWTH change to a better selection proportional to parittion widths  divide big subinterval
      //     i=rand_selection(partition, growth)(0);
        i=randi<int>( distr_param(0, M));
        newpart= newpartitionOPT( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW  ,q,  Pgrams,periodmat,  omegas,  sigma);
        vec psinew=newpart[1];
        vec np=newpart[0];
        
      }
      
    }
  }
  
  mat Weightsnew=GetPhi_part(Z_cur,P_cur,Lambda_prop, newpart[0]);
  mat Weightsold=GetPhi_part(Z_cur,P_cur,Lambda_cur,partition);
  
  double proplhood=Whittle_single_penalized( Pgrams, periodmat,V_cur, Z_cur, newpart[0] , newpart[1], newpart[2],  omegas, sigma );
  double newlhood=0;
  
  accept_Lambda = accept_Lambda + proplhood  -lhood + .01*pow( Lambda_cur,2)-.01*pow(  Lambda_prop,2)  ;  // 
  
  //accept or reject proposed Lambda
  if(accept_Lambda > log(randu())){
    Lambda_new=Lambda_prop;
    newp=as<vec>(newpart[0]);
    newpsi=as<vec>(newpart[1]);
    newBW=as<vec>(newpart[2]);
    newlhood=proplhood;
  }else{
     Lambda_new=Lambda_cur; 
    newlhood=lhood;
    }
  //  }
  
  return List::create(
    _["Lambda_new"] = Lambda_new,
    _["partition"] = newp,
    _["psi"] = newpsi,
    _["BW"] = newBW,
    _["lhood"] = newlhood
    ) ;
}      








//Function to propose new element of V (weights) or Z (atoms) using uniform proposal
// [[Rcpp::export]]
mat NewVZ1(const mat& VZ,const double& epsilon){
  //VZ is matrix of current weights or atoms
  //epsilon is a vector of uniform halfwidths for proposals
  mat prop = VZ +  epsilon*2.0*randu<mat>( size(VZ) ) - ones( size(VZ) )*epsilon;
    //prop.clamp(0, 1);
    mat propB = clamp(prop, 0,     1); 
return (propB );
}



//Function to propose new element of V (weights) or Z (atoms) using uniform proposal
// simplified version
// [[Rcpp::export]]
rowvec NewVZ(const rowvec& VZ,const vec& epsilon){
  
  //VZ is vector of current weights or atoms
  //epsilon is a vector of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
//  double epsi=epsilon(0);
//  mat prop = VZ +  epsi*2.0*randu<mat>( size(VZ) ) - ones( size(VZ) )*epsi;
    rowvec prop=VZ+ epsilon*(2*randu()-1) ;
    return prop;
  
  
  
}


//function to sample stick breaking weights, V, using MH step and uniform proposal
// is going to sample a new row
// [[Rcpp::export]]
List sample_V( const vec& omegas, const vec& psi_cur, const vec& BW_cur,  const cx_mat& Pgrams, mat periodmat, const mat& V_cur, 
              const mat& P_cur,  const mat& Z_cur, int Lambda_cur, const double& epsilon, 
              double alpha_L_cur, double MaxFreq, double sigma, double lhood, vec partition){
  //N is the number of the node to update so the update is univariate and individual
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for time series
  //V_cur is a vector containing the latest iteration BDP stick breaking weights
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column for each subject
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //alpha_L_cur is the latest value of the alpha_L concentration parameter
  //MaxFreq is the maximum frequency to consider for the spectral density
  //DBDP is a boolean to determine if DP atoms are dependent n covariates not implemented 
 // int L=V_cur.n_cols;
 int flag_accept=0;
  int N=V_cur.n_rows;
  mat V_propmat= V_cur;
  mat V_subs= NewVZ1(V_cur, epsilon);
 mat V_newmat=V_cur ;
  double accept_V;  
  //loop through updating each of the L-1 random values of stick breaking weights V
double w0=lhood;
double wN;
  for(int l=0; l<N; l++){
      V_propmat.row(l)=V_subs.row(l);
  wN=Whittle_single_penalized(Pgrams,periodmat, V_propmat, Z_cur , partition   , psi_cur, BW_cur, omegas, sigma );
  accept_V = wN+(alpha_L_cur-1)*sum(log(1-V_propmat.row(l)))  
    - w0-(alpha_L_cur-1)*sum(log(1-V_newmat.row(l))) ;
  //accept or reject new V value using MH step
  if(accept_V > log(randu())){
    V_newmat.row(l)=V_propmat.row(l);
    w0=wN;
    flag_accept=1;
    
  } else{ V_propmat.row(l)= V_newmat.row(l);  flag_accept=0; }
  
  }

  return List::create(
    _["V_new"] = V_newmat,
    _["P_new"] = GetWeightsMult(V_newmat),
    _["lhood"] = w0,
    _["flag_accept"]=flag_accept
  ); 
}


//function to sample atoms, Z, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_Z(  const vec& omegas, const vec&  psi_cur, const vec&  BW_cur, const cx_mat& Pgrams, mat periodmat , const mat& V_cur, 
             const mat& Z_cur, int Lambda_cur, const double& epsilon, double MaxFreq, double sigma, double lhood, vec partition){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a vector containing the latest iteration BDP atoms
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  //lhood id the old value of the whitttle lohlikelihood

int flag_accept=0;
  int N=Z_cur.n_rows ;
  mat Z_newmat=Z_cur;
  mat Z_propmat=Z_cur ;
  mat Z_subs=  NewVZ1(Z_cur, epsilon);
  double accept_Z;  
  double w0=lhood;
  double wN;
 
  for(int l=0; l<N; l++){
    Z_propmat.row(l)=Z_subs.row(l) ; 
    //calculate acceptance probability
    wN=Whittle_single_penalized(Pgrams, periodmat,V_cur,Z_propmat, partition   , psi_cur, BW_cur, omegas, sigma );
    accept_Z = wN- w0 ;
    //accept or reject new V value using MH step
    if(accept_Z > log(randu())){
      Z_newmat.row(l)=Z_propmat.row(l);
      //   V_new_lik=V_prop_lik;
      w0=wN;
      flag_accept=1;
    } else{ Z_propmat.row(l)=Z_newmat.row(l); flag_accept=0; }
  }
  return List::create(
    _["Z_newmat"] = Z_newmat,
    _["lhood"] = w0,
    _["flag_accept"]=flag_accept
  ); 
  }



//function to sample location parameters psi,  given the number of components, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_psi(int N, int L, const vec& psi_cur, const vec& BW_cur,  const vec& omegas, const cx_mat& Pgrams,mat periodmat, const mat& V_cur, const mat& Z_cur,
               int Lambda_cur,  double MaxFreq, vec partition, double sigma , double lhood ){
  //N is the component to generate a new location parameter psi 
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //psi_cur is a vector containing the latest iteration of the location parameters
  //Bwidth is a vector containing the latest iteration of the bandwidth (logmodulus) parameters
  //Lambda_cur is the current number of spectral components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  // lhood is the old value of the whittle loglikelihood
  //compute weights from the P weights of the BDP sampling
  
  double epsilon;
  vec psi_prop(Lambda_cur);
  vec psi_new=psi_cur;
  psi_prop=psi_cur;

  double psi_prop_lik;
  double accept_psi;  

  double psi_new_lik=lhood;
    //propose new location for the l-th subinterval
    epsilon=(partition(N+1)-partition(N))/2.1;
    psi_prop=Newpsij( psi_prop , epsilon, N, partition);
    
  psi_prop_lik = Whittle_single_penalized(Pgrams,periodmat, V_cur,Z_cur,partition, psi_prop, BW_cur, omegas, sigma);  
  //calculate acceptance probability
  accept_psi = psi_prop_lik - lhood;
  //acceptance or reject using MH step
  if(accept_psi > log(randu())){psi_new = psi_prop; psi_new_lik = psi_prop_lik;}
  
  return List::create(
    _["psi_new"] = psi_new,
    _["lhood"] = psi_new_lik
  ) ; 
  
  
}


//function to sample Bandwidth (log modulus) parameters L, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_L(int N, int L, const vec& psi_cur, const vec& BW_cur,  const vec& omegas, const cx_mat& Pgrams, mat periodmat, const mat& V_cur, const mat& Z_cur,
             int Lambda_cur, const double& epsilon, double MaxFreq, double q, double sigma, double lhood, vec partition){
  //N is the component to sample a new bandwidth L
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //psi_cur is a vector containing the latest iteration of the location parameters
  //BW_cur is a vector containing the latest iteration bandwidth (logmodulus) parameters
  //Lambda_cur is the current number of spectral components
  //epsilon is a vector of uniform halfwidths for proposals
  // lhood is the old value of the whittle loglikelihood
  //MaxFreq is the maximum frequency to consider for the spectral density
  //mat Phi=GetPhi_part( Z_cur, P_cur,  Lambda_cur,partition);
  vec BW_prop(Lambda_cur);
  vec BW_new=BW_cur;
  BW_prop=BW_cur;
  int flag_accept=0;
  double BW_new_lik = lhood;
  double BW_prop_lik;
  double accept_BW;  

    //propose new BW
    BW_prop=NewL(BW_prop, epsilon, N, q);
  BW_prop_lik = 0;
  //calculate acceptance probability
  BW_prop_lik = BW_prop_lik +
    Whittle_single_penalized(Pgrams,periodmat,V_cur,Z_cur ,partition, psi_cur, BW_prop, omegas, sigma)   ;
  
  accept_BW = BW_prop_lik - BW_new_lik  + 2*log( BW_prop(N) ) - 2*log( BW_cur(N) ) ;
  //acceptance or reject using MH step
  if(accept_BW > log(randu())){BW_new = BW_prop; BW_new_lik = BW_prop_lik;}
  else{ flag_accept=1; }
  
  return List::create(
    _["BW_new"] = BW_new,
    _["lhood"] = BW_new_lik,
    _["flag_accept"]=flag_accept
  ) ; 
  
}

//Function to propose new variance for the  independent errors 
// [[Rcpp::export]]
double Newsigma(const double& sigma,const double& epsilon){
  //VZ is vector of current weights or atoms
  //epsilon is a vector of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  
  double sigmanew=sigma;
  //propose new value
  double prop=sigma+epsilon*(2*randu()-1);
  //take modulus 1 to ensure it falls between 0 and 1
  if(prop < 0){sigmanew=-prop;}else{sigmanew=prop;}
  return sigmanew;
}

//Function to sample new variance for the  independent errors 
// [[Rcpp::export]]
List sample_sigma( const vec& psi_cur, const vec& BW_cur,  const vec& omegas, const cx_mat& Pgrams, mat periodmat, const mat& V_cur, const mat& Z_cur,
             int Lambda_cur, const double& epsilon, double sigma_cur, double lhood, vec partition ){
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //psi_cur is a vector containing the latest iteration of the location parameters
  //BW_cur is a vector containing the latest iteration bandwidth (logmodulus) parameters
  //Lambda_cur is the current number of spectral components
  //epsilon is a vector of uniform halfwidths for proposals
  //sigma_cur is the current value of the error term variance square root 
  // lhood is the old value of the whittle loglikelihood
   // mat Phi=GetPhi_part( Z_cur, P_cur,  Lambda_cur, partition);
  
  double sigmaprop;
  double sigmaNEW=sigma_cur;
  double new_lik = lhood;
  double prop_lik;
  double accept_BW;  
   //propose new sigma
  sigmaprop=Newsigma(sigma_cur, epsilon );
  prop_lik =  Whittle_single_penalized(Pgrams,periodmat, V_cur,Z_cur,partition , psi_cur, BW_cur, omegas, sigmaprop)   ;
  //calculate acceptance probability
  accept_BW = prop_lik - 2*log( sigmaprop) - new_lik+ 2*log( sigma_cur);
  //acceptance or reject using MH step
  if(accept_BW > log(randu())){
    sigmaNEW = sigmaprop; 
    new_lik=prop_lik;
    }
  
  return List::create(
    _["sigmaNEW"] = sigmaNEW,
    _["lhood"] = new_lik
  ) ; 
  
}



//function to calculate posterior curve samples from posterior variables
// [[Rcpp::export]]
List GetCurvesBDP(int Nsamp, const vec& omegas,  const mat& psi_cur, const mat& BW_cur, const vec& Lambda_samps, 
                  const cube& w){
  //Nsamp is the number of samples
  //omegas is a vector of fourier frequencies
  //L is BDP truncation level
  //Lambda_samps is a vector of samples of the number of Bernstein polynomial components
  //P_samps is a matrix where each row is an L length sample of BDP weights
  //Z_samps is a matrix where each row is an L length sample of BDP atoms
  //MaxFreq is the maximum frequency to consider for the spectral density
  int m =w.n_rows ;
  List CurveEsts(Nsamp);
  mat splinesmatrix ;
  int nOm=omegas.n_elem;
  cube spectralcubes(m,m, nOm,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    cube specmat=spectralmatrix( psi_cur.row(s).t(),BW_cur.row(s).t(), omegas) ;
    spectralcubes=spectralmatrix_nodes( omegas,  w.slice(s), specmat) ;
    CurveEsts[s]=spectralcubes ; 
  }
  return CurveEsts;
}

//%******************* PickAlpha ******************************

// this function computes the posterior distribution  for  the parameter alpha from the Dirichlet Process
// [[Rcpp::export]]
double  getLogPostAlpha(double x, int iStar, double a0, double b0,int n){
  // n is the sample size , for now will be the original sample size of the time series 
  // until new revision for detemrine how it works the histogram maping VS the DFT maping 
  double pi=3.14159265358979323846;
  double alpha = exp(x);
  double h= R::beta(alpha +1, n);
  double logLike = (iStar-1)*log(alpha) +  log(h) + log(alpha+n) ; 
  double logPrior = -.5*log(2*pi)-log(b0)-  .5*(   pow( x - a0,2 )/ pow( b0,2)       )  ;
  double logPost = logLike + logPrior;
  return  logPost ;
}


// This function updates alpha based on the old alpha and number of unique components (slice sampler)
//  
// [[Rcpp::export]]
double  pickAlpha(double oldAlpha, int iStar, double a0, double b0, int n){
  //iStar is the current number of components 
  int  m = 40;
  int  w = 5;
  //   % slice sampling
  double     x = log(oldAlpha+0.001); 
  double    z = getLogPostAlpha(x, iStar, a0, b0,n) +  log(randu()) ;
  //   the exponential ischanged by the log of a a uniform
  double u = randu();
  double   L = x - w*u;
  double   R = L + w;
  double   v = randu();
  int   J = floor(m*v);
  int   K = (m-1) - J;
  
  while(  J>0 && z < getLogPostAlpha(L, iStar, a0, b0, n)){
    L = L - w;
    J = J - 1;}
  
  while( K>0 && z < getLogPostAlpha(R, iStar, a0, b0,n) ){
    R = R+w;
    K = K-1;}
  
  u = randu() ;
  double  newX = L + u*(R-L);
  
  while (z > getLogPostAlpha(newX, iStar, a0, b0,n) ){ 
    if (newX < x){
      L = newX;}
    else{
      R = newX;}
    u = randu();
    newX = L + u*(R-L);
  }
  return  exp(newX);
}  



// Another way to compute the posterior distribution of apha based on a gamma distribution
// [[Rcpp::export]]
double  alphamixgampost( double oldeta, int oldlambda , double a0, double b0, int n){
double newalpha;  
double x= (a0 +oldlambda-1)/(static_cast<double>(n)*(b0-log(oldeta)  )  );
double pieta=x/(1+x) ;  
double u =randu();

if( pieta<=u ){
  newalpha= R::rgamma( a0+oldlambda,  pow(b0-log(oldeta),-1)    )  ;
}else{
  newalpha= R::rgamma(a0+oldlambda-1, pow(b0-log(oldeta),-1)   );
     }
return(newalpha)   ; 
}



//MCMC Sampling Algorithm for the MBMARD model
// [[Rcpp::export]]
List MBMARD(double  Tsize ,const vec& omegas,const cx_mat& Pgram,mat periodmat, int Nsamp, int L, const double& epsilonV, const double& epsilonZ,  const double& epsilon_BW, mat Zinit, mat Vinit, vec psi_init, vec BWinit,  
                 double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                 int Lambda_init=30, string Lambda_prior="flat", int Lambda_max=300, 
                 string Lambda_prop_dist="poisson", bool alpha_const=true, 
                 double alpha_init=1, double a_alpha_L=1, double b_alpha_L=1  , double q=1, double sigma=1, double epsilon_sigma=.1    ){
  //X is a vector containing the time series to be analyzed
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //epsilon_BW is a vector of uniform halfwidths for Bandwidth proposals length same as lambda_init
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
  //either "expon", "poisson", or "flat"
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps
  //either "poisson" or "up1down1". If mispecifiec or "const" then lambda will be constant 
  //alpha_const is a boolean to determine if the concetration parameter is set to constant
  //alpha_init is the initial value for concentration parameter and the value for all samples if constant
  //a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha is to be sampled
  //q is the limit for the uniform distribution of the bandwidth parameters


  vec sigmavec(Nsamp) ;
  sigmavec(0)=sigma ;
  //number of nodes
  int m= Pgram.n_rows;
  //initialize concentration parameter alpha
  vec alpha_L(Nsamp);
  alpha_L(0)=alpha_init;

  //initialize BDP (from Choudhuri)
  vec Lambda(Nsamp);
  Lambda(0)=Lambda_init;
  cube w(m ,Lambda_max,Nsamp, fill::zeros);
  mat V(m, L);
  V=Vinit    ;
  mat P(m, L);
  P=GetWeightsMult(V);
  List NewVP(3);
  List Ztemp(2);
  List sigmatemp(2);
  List bwtemp(2);
  List psitemp(2);
  List  Birth_death(5);
  mat Z(m, L);
  Z =Zinit;
  double  templhood ;
  //random start for each parameter in order to evaluate the G-R psrf 
  double eta;
  mat psi(Nsamp, Lambda_max, fill::zeros);
  mat BW(Nsamp, Lambda_max);
  BW.row(0).subvec(0,Lambda_init-1)=BWinit.t();
  psi.row(0).subvec(0,Lambda_init-1)=psi_init.t(); 
  mat partition(Nsamp, Lambda_max+1, fill::zeros);
  vec partition_init= zeros<vec>(Lambda_init+1)  ; 
  // inititalization of the parition 
  partition_init(Lambda_init)=0.5 ;

  for(int i=1 ; i<Lambda_init ;i++){
    partition_init(i)= (psi_init(i-1) + psi_init(i) )/2.0    ;
  }  
  
  partition.row(0).subvec(0,Lambda_init)=partition_init.t();
  
  w.subcube(0,0,0,m-1,Lambda_init-1,0)=GetPhi_part(Z, P, Lambda_init,partition_init);

  vec loglik(Nsamp);
  loglik(0)=Whittle_single_penalized(Pgram, periodmat,Z,V,partition_init  , psi.row(0).subvec(0,Lambda_init-1).t() , BW.row(0).subvec(0,Lambda_init-1).t(), omegas, sigmavec(0)  );
  
  
// iterations for the rest of the algorithm  
  
   for(int iter=1; iter<Nsamp; iter++){

    //Update number of Bernstein polynomials, Lambda
    Birth_death     = sample_Lambda(1, L, omegas, psi.row(iter-1).subvec(0,Lambda(iter-1)-1).t() , BW.row(iter-1).subvec(0,Lambda(iter-1)-1).t(), Pgram, periodmat, P, V ,  Z, epsilon_BW,
                                          Lambda(iter-1), Lambda_prior, Lambda_max, MaxFreq,  partition.row(iter-1).subvec(0,Lambda(iter-1)).t(), q, sigmavec(iter-1), loglik(iter -1));

    Lambda(iter)=as<int>( Birth_death[0]);
    partition.row(iter).subvec(0,Lambda(iter))   =as<vec>(Birth_death[1]).t();

  // to maka an intermediate row and update the psi and bwand update the intermediate,
      templhood =as<double>(Birth_death[4]);
      psi.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(Birth_death[2]).t() ;
      BW.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(Birth_death[3]).t();

    for(int l=0; l<Lambda(iter); l++){
      // Update psi
    //list sample psi   with previous temp hood
    // then psi row
     psitemp = sample_psi(l,  L,  psi.row(iter).subvec(0,Lambda(iter)-1).t(), BW.row(iter).subvec(0,Lambda(iter)-1).t(), omegas, Pgram,periodmat, V, Z,
            Lambda(iter)     ,  MaxFreq,partition.row(iter).subvec(0,Lambda(iter)).t() , sigmavec(iter-1), templhood   );
      psi.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(psitemp[0]).t();
      //list sample bw
      templhood=as<double>(psitemp[1] );
      // then bw row
      //Update Bandwidths
       bwtemp= sample_L(l,  L,   psi.row(iter).subvec(0,Lambda(iter)-1).t() , BW.row(iter).subvec(0,Lambda(iter)-1).t() , omegas, Pgram,periodmat, V, Z,
             Lambda(iter), epsilon_BW,  MaxFreq, q, sigmavec(iter-1),templhood, partition.row(iter).subvec(0,Lambda(iter)).t() );
      BW.row(iter).subvec(0,Lambda(iter)-1)=as<vec>( bwtemp[0]).t();
      templhood=as<double>(bwtemp[1]);
          }
// updates of sigma
     sigmatemp=sample_sigma(psi.row(iter).subvec(0,Lambda(iter)-1).t() , BW.row(iter).subvec(0,Lambda(iter)-1).t(),  omegas,  Pgram,periodmat, V, Z ,
                                Lambda(iter), epsilon_sigma, sigmavec(iter-1),templhood ,partition.row(iter).subvec(0,Lambda(iter)).t() ) ;
    sigmavec(iter)= as<double>(sigmatemp[0]) ;

    // after sampling the number of components the weights are obtained again using  the same atoms and bdp weigths
    //the partition is proposed selecting one subinterval at random
    // then a new point is drawn randomly on the interval selected
    //  a new location location is drawn in he subpartition where there is no location parameters  psi_jump
    // a new bandwidth is proposed uniformly
    // in the death process parameters are substituted by one value drawn randomly according to the transition probabilities
    //in location the mid point of the two parameters is select and from it a new location is proposed in
    // the selected subinterval

    //Update weights V and P

    NewVP = sample_V( omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(), BW.row(iter).subvec(0,Lambda(iter)-1).t(), Pgram,periodmat, V, P,
                     Z, Lambda(iter), epsilonV, alpha_L(iter-1), MaxFreq, sigmavec(iter),as<double>( sigmatemp[1]) ,partition.row(iter).subvec(0,Lambda(iter)).t() );
    V = as<mat>( NewVP("V_new")) ;
    P = as<mat>( NewVP("P_new") );

    //Update Z

     Ztemp= sample_Z( omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(),BW.row(iter).subvec(0,Lambda(iter)-1).t(),  Pgram,periodmat, V, Z,
                          Lambda(iter), epsilonZ, MaxFreq, sigmavec(iter),as<double>(NewVP[2]),partition.row(iter).subvec(0,Lambda(iter)).t() );
    Z = as<mat>( Ztemp[0] );

    //Dirichlet process Update alpha_L via Gibbs step
    // modofication gibss build from gamma prior west 1995
    eta= R::rbeta( alpha_L(iter-1)+1, Tsize);

    if(alpha_const==false){
      //   alpha_L(iter)=my_rgamma(a_alpha_L+L-1, b_alpha_L - log(P(iter,L-1))    )[0] ;
      //this update is based on Ishwaran and James 2002 ,
      // based on the truncated DP approximation
      // alpha | V /sim GA( L +a -1, b- log pm)  pm is the last weight of the truncated atoms
      // update below is based on slice sapling and the gibbs sampler from west and escobar  (1995 ), method from prof shahbaba
    //  alpha_L(iter)=pickAlpha(alpha_L(iter-1), Lambda(iter) ,  a_alpha_L  , b_alpha_L, Tsize  );
    alpha_L(iter)=alphamixgampost( eta, Lambda(iter), a_alpha_L, b_alpha_L, Tsize);

    }
    if(alpha_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }

    //optional sentence to send a message each Sup iterations
   // if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}

    // retrieve the weights
    w.subcube(0,0,iter,m-1,Lambda(iter)-1,iter)=GetPhi_part(Z, P, Lambda(iter) ,partition.row(iter).subvec(0,Lambda(iter)).t() );

    loglik(iter)= as<double>(Ztemp[1]);
   }

  // optional return the spectral curves of each iterations, outcomes could be Large
  //List Curves=GetCurvesBDP(Nsamp, omegas,psi,BW,   Lambda, w);

  // return only effective samples, select the parameter you would like to monitor
  return List::create(
    _["Lambda"] = Lambda,
    _["weights"] = w,
    _["psi"] = psi,
    _["BW"] = BW,
  //  _["CurveEsts"] = Curves,
    _["V"] = V   ,
    _["Z"] = Z,
    _["omegas"] = omegas,
    _["Pgram"] = Pgram ,
   _["alpha"] = alpha_L,
   _["Whittle"] = loglik,
   _["Tsize"] = Tsize,
//   _["TruncL"] = L ,
   _["partition"] = partition.row(Nsamp-1) ,
   _["sigma"] = sigmavec
//   _["epsilonZ"] = epsilonZ,
//   _["epsilonBW"] = epsilon_BW,
//   _["SamplingRate"] = SamplingRate,
//   _["MaxFreq"] = MaxFreq,
//   _["Lambda_max"] = Lambda_max,
//   _["alpha_const"] = alpha_const,
//   _["a_alpha_L"] = a_alpha_L,
//   _["b_alpha_L"] = b_alpha_L,
//   _["q"] = q
  );

  

  
  
}

//MCMC Sampling from a started MBMARD chain
// [[Rcpp::export]]
List MBMARD_continuation(  cx_mat Pgram,mat periodmat , 
                                  mat Vold,
                                  mat Zold ,
                                  double alphaold,
                                  int Lamold ,
                                  mat wold ,
                                  vec psiold  ,
                                  vec BWold  ,
                                  vec omegas,
                                  double whittleold ,
                                  vec partitionold  ,
                                  double sigmaold,
                                  int Tsize , int Nsamp,  const double& epsilonV, const double& epsilonZ,  const double& epsilon_BW,  
                 double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                  string Lambda_prior="flat",  
                 string Lambda_prop_dist="poisson", bool alpha_const=true, 
                  double a_alpha_L=1, double b_alpha_L=1  , double q=1 , double epsilon_sigma =.01  ){
  //X is a vector containing the time series to be analyzed
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //epsilon_BW is a vector of uniform halfwidths for Bandwidth proposals length same as lambda_init
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
  //either "expon", "poisson", or "flat"
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps
  //either "poisson" or "up1down1". If mispecifiec or "const" then lambda will be constant 
  //alpha_const is a boolean to determine if the concetration parameter is set to constant
  //alpha_init is the initial value for concentration parameter and the value for all samples if constant
  //a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha is to be sampled
  //q is the limit for the uniform distribution of the bandwidth parameters


vec sigmavec(Nsamp) ;
  sigmavec(0)=sigmaold ;
  int L=Zold.n_cols ;  // size of truncation
  //number of nodes
  int m= Pgram.n_rows;
  //initialize concentration parameter alpha
  vec alpha_L(Nsamp);
  alpha_L(0)= alphaold ;
  vec Lambda(Nsamp);
  Lambda(0)=Lamold;
  cube w(m ,wold.n_cols,Nsamp, fill::zeros);
  mat V=Vold;
  mat P(m, L);
  P=GetWeightsMult(V);
  List NewVP(3);
  List Ztemp(2);
  List sigmatemp(2);
  List bwtemp(2);
  List psitemp(2);
  List  Birth_death(5);
  
  mat Z=Zold;
   //random start for each parameter in order to evaluate the G-R psrf 
  double eta;
  mat psi(Nsamp, psiold.size(), fill::zeros);
  mat BW(Nsamp, BWold.size());
  BW.fill(   randu()*q   );
  BW.row(0)=BWold.t();
  mat partition(Nsamp, partitionold.size(), fill::zeros);
    psi.row(0)=psiold.t(); 
  partition.row(0)=partitionold.t() ;
  double templhood;
  w.slice(0)=wold;
  vec loglik(Nsamp);
  loglik(0)=whittleold;
  for(int iter=1; iter<Nsamp; iter++){
    //Update number of Bernstein polynomials, Lambda
    List  Birth_death     = sample_Lambda(1, L, omegas, psi.row(iter-1).subvec(0,Lambda(iter-1)-1).t() , BW.row(iter-1).subvec(0,Lambda(iter-1)-1).t(), Pgram,periodmat, P, V ,  Z, epsilon_BW,
                                          Lambda(iter-1), Lambda_prior,  wold.n_cols, MaxFreq,  partition.row(iter-1).subvec(0,Lambda(iter-1)).t(), q, sigmavec(iter-1),  loglik(iter-1) );
    Lambda(iter)=as<int>( Birth_death[0]);
    
    partition.row(iter).subvec(0,Lambda(iter))   =as<vec>(Birth_death[1]).t();
    templhood =as<double>(Birth_death[4]);
    
    psi.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(Birth_death[2]).t() ;
    BW.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(Birth_death[3]).t();
    
    for(int l=0; l<Lambda(iter); l++){
      // Update psi
      //list sample psi   with previous temp hood
      //temp lhood
      // then psi row  
      psitemp = sample_psi(l,  L,  psi.row(iter).subvec(0,Lambda(iter)-1).t(), BW.row(iter).subvec(0,Lambda(iter)-1).t(), omegas, Pgram,periodmat, V, Z,
                           Lambda(iter)     ,  MaxFreq,partition.row(iter).subvec(0,Lambda(iter)).t() , sigmavec(iter-1), templhood   );
      psi.row(iter).subvec(0,Lambda(iter)-1)=as<vec>(psitemp[0]).t(); 
      //list sample bw
      templhood=as<double>(psitemp[1] );
      // then bw row  
      //Update Bandwidths
      
      bwtemp= sample_L(l,  L,   psi.row(iter).subvec(0,Lambda(iter)-1).t() , BW.row(iter).subvec(0,Lambda(iter)-1).t() , omegas, Pgram,periodmat, V, Z,
                       Lambda(iter), epsilon_BW,  MaxFreq, q, sigmavec(iter-1),templhood,partition.row(iter).subvec(0,Lambda(iter)).t()  );
      BW.row(iter).subvec(0,Lambda(iter)-1)=as<vec>( bwtemp[0]).t();
      templhood=as<double>(bwtemp[1]);
      
    }
    
    
  // sample sigma   
    sigmatemp=sample_sigma(psi.row(iter).subvec(0,Lambda(iter)-1).t() , BW.row(iter).subvec(0,Lambda(iter)-1).t(),  omegas,  Pgram,periodmat, V, Z,
                           Lambda(iter), epsilon_sigma, sigmavec(iter-1),templhood ,partition.row(iter).subvec(0,Lambda(iter)).t() ) ;
    
    sigmavec(iter)= as<double>(sigmatemp[0]) ;     
    
    //Update weights V and P
    
    NewVP = sample_V( omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(), BW.row(iter).subvec(0,Lambda(iter)-1).t(), Pgram,periodmat, V, P,
                      Z, Lambda(iter), epsilonV, alpha_L(iter-1), MaxFreq, sigmavec(iter),as<double>( sigmatemp[1]),partition.row(iter).subvec(0,Lambda(iter)).t()   );
    V = as<mat>( NewVP("V_new")) ;
    P= as<mat>( NewVP("P_new") );
    
    //Update Z
    
    Ztemp= sample_Z( omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(),BW.row(iter).subvec(0,Lambda(iter)-1).t(),  Pgram,periodmat, V, Z,
                     Lambda(iter), epsilonZ, MaxFreq, sigmavec(iter),as<double>(NewVP[2]),partition.row(iter).subvec(0,Lambda(iter)).t()  );
    
    Z = as<mat>( Ztemp[0] );
    
    //Dirichlet process Update alpha_L via Gibbs step
    
    // modofication gibss build from gamma prior west 1995
    eta= R::rbeta( alpha_L(iter-1)+1, Tsize);

    if(alpha_const==false){
      alpha_L(iter)=alphamixgampost( eta, Lambda(iter), a_alpha_L, b_alpha_L, Tsize);
    }
    if(alpha_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }
    
    
    //optional message
    // if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}
    
    // retrieve the weights
    w.subcube(0,0,iter,m-1,Lambda(iter)-1,iter)=GetPhi_part(Z, P, Lambda(iter),partition.row(iter).subvec(0,Lambda(iter)).t() );
    loglik(iter)= as<double>(Ztemp[1]);
  }
// optional computation if the spectral matrices
    //List Curves=GetCurvesBDP(Nsamp, omegas,psi,BW,   Lambda, w);
  
  // return only effective samples
  return List::create(
    _["Lambda"] = Lambda,
    _["weights"] = w,
    _["psi"] = psi,
    _["BW"] = BW,
  //  _["CurveEsts"] = Curves,
    _["V"] = V   ,
    _["Z"] = Z,
    _["omegas"] = omegas,
    _["Pgram"] = Pgram ,  
    _["alpha"] = alpha_L, 
    _["Whittle"] = loglik, 
    _["Tsize"] = Tsize,
    //   _["TruncL"] = L ,
    _["partition"] = partition.row(Nsamp-1) ,
    _["sigma"] = sigmavec
    //   _["epsilonZ"] = epsilonZ,
    //   _["epsilonBW"] = epsilon_BW,
    //   _["SamplingRate"] = SamplingRate,
    //   _["MaxFreq"] = MaxFreq,
    //   _["Lambda_max"] = Lambda_max,
    //   _["alpha_const"] = alpha_const,
    //   _["a_alpha_L"] = a_alpha_L,
    //   _["b_alpha_L"] = b_alpha_L,
    //   _["q"] = q
  );
  
}










