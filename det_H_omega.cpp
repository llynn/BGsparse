// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <Rcpp/Benchmark/Timer.h>

// [[Rcpp::export]]
double det_H_omega(arma::mat invOmega, IntegerVector gam, int nu, int p, IntegerMatrix indi, int pC2){
  arma::mat H(nu, nu);
  Timer tt;

  IntegerMatrix indi_sub(nu,2);
  int count = 0; int t = 0;
  for(int s = 0; s < (pC2+p); s ++){
    if(indi(s,0) == indi(s,1)){
      indi_sub(count,_) = indi(s,_);
      count +=1;
    }else{
      if(gam[t] == 1){
        indi_sub(count,_) = indi(s,_);
        count+=1;
      }
      t +=1;
    }
  }
  
  tt.step("indi_sub") ;

  int xi; int xj; int xl; int xm;
  for(int i = 0; i < nu; i++){
    xi = indi_sub(i,0); xj = indi_sub(i,1);
    //diagnoal
    if(xi == xj){
      H(i,i) = -invOmega(xi,xj)*invOmega(xi,xj);
    }else{
      H(i,i) = -2*invOmega(xj,xi)*invOmega(xj, xi) - 2*invOmega(xi, xi)*invOmega(xj, xj);
    }
    for (int j = (i+1); j < nu; j++){
        xl = indi_sub(j,0); xm = indi_sub(j,1);
        if((xi != xj)&(xl != xm)) {
          H(i,j) = -invOmega(xj,xl)*invOmega(xm, xi) - invOmega(xi, xl)*invOmega(xm, xj) - 
                 invOmega(xj,xm)*invOmega(xl, xi) - invOmega(xi, xm)*invOmega(xl, xj);
          
        }else{
          if((xl == xm) & (xi == xj)){
            H(i,j) = -invOmega(xi,xl)*invOmega(xl, xi);
            
          }else if(xi == xj){
            H(i,j) = -invOmega(xj, xl)*invOmega(xm,xi)-invOmega(xi, xm)*invOmega(xl,xj);
            
          }else{
            H(i,j) = -invOmega(xj, xl)*invOmega(xl,xi)-invOmega(xi, xl)*invOmega(xl,xj);
          }
            
        }
        H(j,i) = H(i,j);
          
    }
  }
  tt.step("get H") ;
  //Rcout << H<<std::endl;
  
  //arma::vec eigval;
  //arma:: mat eigvec;
  //eig_sym(eigval, eigvec, H);
  //Rcout << eigval <<std::endl;
  
    //Rcout << U<<std::endl;

  double output = det(H);
  tt.step("det H") ;
  
  arma::mat res(REAL(tt),3,1,false,true);
  Rcout << res<<std::endl;
  Rcout << H ;
  return output;
}


// [[Rcpp::export]]
double det_H_omega_slow(arma::mat invOmega, IntegerVector gam, int nu, int p, IntegerMatrix indi, int pC2){
  // arma::mat invOmega(REAL(x), x.nrow(), x.ncol(), false, true);
  arma::mat xE1(p, p);
  arma::mat xE2(p, p);
  arma::mat H(nu, nu);
  Timer tt;
  //Rcout<<"nu = "<<nu << " p = "<<p << " pC2" << pC2 <<std::endl;
//  for(int s = 0; s < pC2; s++){
//    Rcout << "gam[s] = " << gam[s] << " s = " << s << std::endl;
//  }
  IntegerMatrix indi_sub(nu,2);
  int count = 0; int t = 0;
  for(int s = 0; s < (pC2+p); s ++){
    if(indi(s,0) == indi(s,1)){
      indi_sub(count,_) = indi(s,_);
      count +=1;
    }else{
      if(gam[t] == 1){
        indi_sub(count,_) = indi(s,_);
        count+=1;
      }
      t +=1;
    }
  }
  
  tt.step("indi_sub") ;

//  
//  for(int k = 0; k < nu; k++){
//    Rcout<<"k =  "<< k << " indi_sub = "<< indi_sub(k,0) << " , " << indi_sub(k,1) << "indi(k,0) = "<<indi(k,0) <<std::endl;
//  }
  
  for(int i = 0; i < nu; i++){
    for (int j = i; j < nu; j++){
        xE1.zeros(); xE2.zeros();
        xE1(indi_sub(i,0),indi_sub(i,1)) = 1; xE1(indi_sub(i,1),indi_sub(i,0)) = 1;
        xE2(indi_sub(j,0),indi_sub(j,1)) = 1; xE2(indi_sub(j,1),indi_sub(j,0)) = 1;
//        if(i ==0 & j == 1){
//          Rcout << xE1 << std::endl;
//          Rcout << xE2 << std::endl;
//        }
//        
        H(i,j) = -trace(invOmega * xE1 * invOmega * xE2);
        H(j,i) = H(i,j);
        //xE1(indi_sub(i,0),indi_sub(i,1)) = 0; xE1(indi_sub(i,1),indi_sub(i,0)) = 0;
        //xE2(indi_sub(j,0),indi_sub(j,1)) = 0; xE2(indi_sub(j,1),indi_sub(j,0)) = 0;
        //Rcout<<xE1<<std::endl;
    }
  }
  tt.step("get H") ;
  //Rcout << H<<std::endl;
  
  //arma::vec eigval;
  //arma:: mat eigvec;
  //eig_sym(eigval, eigvec, H);
  //Rcout << eigval <<std::endl;
  
    //Rcout << U<<std::endl;

  double output = det(H);
  tt.step("det H") ;
  
  arma::mat res(REAL(tt),3,1,false,true);
  Rcout << res<<std::endl;
  Rcout << H ;
  return output;
}
