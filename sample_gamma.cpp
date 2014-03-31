// sample gamma, also update q 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <Rcpp/Benchmark/Timer.h>

//double trace(NumericMatrix m) {  
//  // compute trace of a square matrix
//  int p = m.ncol();
//  double tr = 0.;
//  for (int i = 0; i < p; i++){
//    tr += m(i,i);
//  }
//  return tr;
//}

double max(double m, double n) {  
  // return max value
  double output;
  if (m > n){
    output = m;
  }else{
    output = n;
  }
  return output;
}

double det_H_omega(arma::mat invOmega, IntegerVector gam, int nu, int p, IntegerMatrix indi, int pC2){
 arma::mat H(nu, nu);
  //Timer tt;
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
  
  //tt.step("indi_sub") ;

  int xi; int xj; int xl; int xm;
  for(int i = 0; i < nu; i++){
    xi = indi_sub(i,0); xj = indi_sub(i,1);
    //diagonal
    if(xi == xj){
      H(i,i) = -invOmega(xi,xj)*invOmega(xi,xj);
    }else{
      H(i,i) = -2*invOmega(xj,xi)*invOmega(xj, xi) - 2*invOmega(xi, xi)*invOmega(xj, xj);
    }
    //off-diagonal
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
  //tt.step("get H") ;

  double output = det(H);
  //tt.step("det H") ;
  
  //arma::mat res(REAL(tt), 3, 1, false, true);
 // Rcout << res<<std::endl;
  //Rcout << H ;
  return output;
}

// [[Rcpp::export]]
SEXP sample_gamma(NumericMatrix xOmega_star_inv, NumericMatrix xS, double xlambda, int xn, IntegerVector xgamma,
                  NumericVector xome_upper, double p_gamma, double pg1, double pg2, 
                  int r_bar, double q, double a_q, double b_q, IntegerMatrix indi) {   
    int pC2 = xome_upper.size();
    int p = xS.ncol();
    int count = 0; // for the length of freemove
    int accept = 0;
    IntegerVector gamma_tt = clone(xgamma);
//    arma::mat Omega_star(REAL(xOmega_star),xOmega_star.nrow(), xOmega_star.ncol(), true, true); 
//    arma::mat U_inv = inv(chol(Omega_star));
//
//    arma::mat Omega_star_inv = U_inv*U_inv.t();
//    timer.step("find inverse") ;
//    NumericVector res(timer);
//    Rcout << res[0] << std::endl;
    
    arma::mat Omega_star_inv(REAL(xOmega_star_inv),p, p, true, true); 
    for(int i = 0; i < pC2; i++){
      if (xome_upper[i] != 0) {
        count +=1;
      }
    }
    IntegerVector fixed_position(pC2-count); int flag1 = 0;
    IntegerVector freemove(count); int flag2 = 0;
    for(int i=0; i<pC2; i++){
      if(xome_upper[i] == 0){
        fixed_position[flag1] = i;
        flag1 +=1;
      }else{
        freemove[flag2] = i;
        flag2 +=1;
      }
    }
    
    //to propose new gamma by randomly select a few positions in freemove, and change the value by 1-gamma
    if(count > 0){

      if(count == 1){
        gamma_tt[freemove[0]] = 1-xgamma[freemove[0]];
      }else{
         Rcpp::NumericVector m(1); 
        do{
          gamma_tt = clone(xgamma);
          if (Rcpp::as<double>(Rcpp::runif(1)) <= p_gamma) {
            m = Rcpp::rbinom(1,count,pg1); 
          } else {
            m = Rcpp::rbinom(1,count,pg2); 
          }
          
          if(m[0] == 0){m[0] +=1;}

          /* random sample */
          int vec = 0; int arr = 0;
          for (int k = 0; k < m[0]; k++){
            NumericVector tmp = runif(1,k,count);
            vec = (int) tmp[0];
            arr = freemove[k];
            freemove[k] = freemove[vec];
            freemove[vec] = arr;
            gamma_tt[freemove[k]] = 1-gamma_tt[freemove[k]];
          }
        }while(sum(gamma_tt)>r_bar);
        
      }
      
      
      // comupte log posterior distribution
      double log_deno = 0; double log_nume = 0;
      double lgamma_tt = sum(gamma_tt); double lgamma = sum(xgamma);
      //Rcout<<"lgamma_tt = "<<lgamma_tt <<"lgamma = "<<lgamma<<std::endl;

      log_nume += lgamma_tt*log(q) + (pC2-lgamma_tt)*log(1-q) + (p+lgamma_tt)*log(xlambda/2) + 0.5*(lgamma_tt + p)*(log(2*M_PI) - log(xn/2));
      log_deno += lgamma*log(q) + (pC2-lgamma)*log(1-q) + (p+lgamma)*log(xlambda/2) + 0.5*(lgamma + p)*(log(2*M_PI) - log(xn/2));

      for( int i=0; i < pC2; i++){
        if(gamma_tt[i] ==1){log_nume -= xlambda*abs(xome_upper[i]);}
        if(xgamma[i] == 1) {log_deno -= xlambda*abs(xome_upper[i]);}
      }
      
      //NumericMatrix H = det_H_omega(Omega_star_inv, gamma_tt, (p+lgamma_tt), p, indi, pC2);
      
      double temp = det_H_omega(Omega_star_inv, gamma_tt, (p+lgamma_tt), p, indi, pC2);
      log_nume -= 0.5*log(max(temp, 1.0e-20));

      //Rcout << "determinant = "<<temp<<std::endl;
      
      temp = det_H_omega(Omega_star_inv, xgamma, (p+lgamma), p, indi, pC2);
      
      log_deno -= 0.5*log(max(temp, 1.0e-20));
      //Rcout << "determinant = "<<temp<<std::endl;
      //Rcout<<"4. log_nume = " << log_nume <<"log_deno = "<<log_deno<<std::endl;
      
       if (log(as<double>(runif(1)) ) <= (log_nume - log_deno)) {
         accept = 1;
       }else{
         gamma_tt = xgamma;
       }
    }
    
    // to update q, restrict q to be less than 0.5
    NumericVector qq = rbeta(1,sum(gamma_tt)+a_q, pC2-sum(gamma_tt)+b_q);
    if(qq[0] < 0.5){
      q = qq[0];   
    }
    
    //return wrap(as<double>(R_det(xOmega_star)));
    return List::create(Named("accept") = accept, Named("gamma_tt") = gamma_tt, Named("q") = q);

}