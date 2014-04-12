// sample gamma, also update q 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <Rcpp/Benchmark/Timer.h>

int sample_index(int n, NumericVector weights) {  
  // normalize weights from log weights
  double m = max(weights);
  double logsum=0.;
  for(int ii = 0; ii < n; ii++){
    logsum += exp(weights[ii] - m);
  }
  logsum = m + log(logsum);
  
  weights = exp(weights-logsum);
  //arma::vec wei = as<arma::vec>(weights);
  
  bool flag = false;
  int index = 0;
  double lower = 0.0;
  double upper = weights[0];
  double randU = as<double>(runif(1));
  
  // find out which interval the sampled number lies within
  do {
    if ((randU >= lower) & (randU < upper)) {
      flag = true;
    } else {
      index +=1;
      lower = upper;
      upper = upper + weights[index]; 
    }
  } while (!flag);
  return index; // start from 0
}


double det_H_omega(arma::mat invOmega, IntegerVector gam, int nu, int p, IntegerMatrix indi, int pC2){
  arma::mat H(nu, nu);
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

  int xi; int xj; int xl; int xm;
  for(int i = 0; i < nu; i++){
    xi = indi_sub(i,0); xj = indi_sub(i,1);
    //diagonal
    if(xi == xj){
      H(i,i) = invOmega(xi,xj)*invOmega(xi,xj);
    }else{
      H(i,i) = 2*invOmega(xj,xi)*invOmega(xj, xi) + 2*invOmega(xi, xi)*invOmega(xj, xj);
    }
    //off-diagonal
    for (int j = (i+1); j < nu; j++){
        xl = indi_sub(j,0); xm = indi_sub(j,1);
        if((xi != xj)&(xl != xm)) {
          H(i,j) = invOmega(xj,xl)*invOmega(xm, xi) + invOmega(xi, xl)*invOmega(xm, xj) + 
                 invOmega(xj,xm)*invOmega(xl, xi) + invOmega(xi, xm)*invOmega(xl, xj);
        }else{
          if((xl == xm) & (xi == xj)){
            H(i,j) = invOmega(xi,xl)*invOmega(xl, xi);
            
          }else if(xi == xj){
            H(i,j) = invOmega(xj, xl)*invOmega(xm,xi)+invOmega(xi, xm)*invOmega(xl,xj);
            
          }else{
            H(i,j) = invOmega(xj, xl)*invOmega(xl,xi)+invOmega(xi, xl)*invOmega(xl,xj);
          }
            
        }
        H(j,i) = H(i,j);
          
    }
  }

  double output = det(H);
  return output;
}

// [[Rcpp::export]]
SEXP sample_gamma_exact(NumericMatrix xOmega_star_inv, NumericMatrix xOmega_star, NumericMatrix SS, double xlambda, int xn, 
                  NumericVector xome_upper, double q, double a_q, double b_q, IntegerMatrix indi, IntegerMatrix d) {   
    int pC2 = xome_upper.size();
    int p = SS.ncol();
    IntegerVector gamma_tt(pC2);
    IntegerVector gamma_max(pC2);

    arma::mat Omega_star_inv(REAL(xOmega_star_inv),p, p, false, true); 
    arma::mat Omega_star(REAL(xOmega_star),p, p, false, true); 
    arma::mat xS(REAL(SS),p, p, false, true); 
  
    int drow = d.nrow();
    NumericVector weights(drow); //log posteri
    
    double constant =  xn*(log(det(Omega_star)) - trace(xS*Omega_star))/2 - xlambda*trace(Omega_star)/2; // eqn 4.2 for Omega*
    
    //to calculate the distribution of gamma for all possible combination
    if(drow >0){
      for(int ii = 0; ii < drow; ii ++){
        int lgamma = sum(d(ii,_));
        weights[ii] =  lgamma*log(q) + (pC2-lgamma)*log(1-q) + (p+lgamma)*log(xlambda/2) + 0.5*(lgamma + p)*(log(2*M_PI) - log(xn/2)); //
        for( int i=0; i < pC2; i++){
          if(d(ii,i) ==1){weights[ii] -= xlambda*abs(xome_upper[i]);}
        }
        double temp = det_H_omega(Omega_star_inv, d(ii,_), (p+lgamma), p, indi, pC2);
        weights[ii] -= 0.5*log(temp);
        weights[ii] += constant;
      }
      
      //arma::vec wei = as<arma::vec>(weights);
      //Rcout<<wei;
      int index = sample_index(drow, weights);
      gamma_tt = d(index,_);
      
      int maxm = 0;
      for(int ii = 1; ii < drow; ii++){
        if(weights[ii] > weights[maxm]) {maxm = ii;}
      }
      gamma_max = d(maxm,_);
    }
    
    // to update q, restrict q to be less than 0.5
    NumericVector qq = rbeta(1,sum(gamma_tt)+a_q, pC2-sum(gamma_tt)+b_q);
    q = qq[0];
//    if(qq[0] < 0.5){
//      q = qq[0];   
//    }
    
    //return wrap(as<double>(R_det(xOmega_star)));
    return List::create( Named("gamma_tt") = gamma_tt, Named("q") = q, Named("gamma_max") = gamma_max, Named("weights") = weights);

}