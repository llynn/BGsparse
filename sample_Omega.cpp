// sample Omega columnwise, together with lambda

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP sample_Omega(NumericMatrix xOmega_star, NumericMatrix xS, IntegerMatrix ind_noi_all, double xlambda, int xn, 
                  IntegerMatrix gamma_M, NumericVector xome_upper, NumericMatrix tau, NumericMatrix xSig, IntegerMatrix indi,
                  NumericMatrix Omega, double a_post, double b_lambda) {
   int p = xS.ncol();
   IntegerVector ind_noi(p-1);
   arma::mat Sig(REAL(xSig), xSig.nrow(), xSig.ncol(), false, true);
   arma::mat xOmega(REAL(Omega), Omega.nrow(), Omega.ncol(), false, true);
    
   double b_post;
   NumericVector gam;
   for (int i = 0; i < p; i ++){
     gam = rgamma(1, xn/2.0+1, 2.0/(xS(i,i)+xlambda));
     xOmega(i,i) = gam[0];
       
     ind_noi = ind_noi_all(_,i);
     int count = 0;
     for(int s = 0; s < p-1; s ++){
       if(gamma_M(ind_noi[s],i) == 1){count +=1;}
     }
     if(count > 0){    
     arma::mat Sig11(count, count);
     arma::mat Sig12(count,1);
     arma::mat invC11(count,count);
     arma::mat Ci(count, count);
     arma::mat Ci_chol(count, count);
     arma::mat diagonal(count,count);
     arma::colvec tau_temp(count,1);
     arma::mat S12(count,1);
     
     count = 0;
     for(int s = 0; s < p-1; s++){
       if(gamma_M(ind_noi[s],i) == 1){
         tau_temp[count] = tau(ind_noi[s],i);
         Sig12[count] = Sig(ind_noi[s],i);
         S12[count] = xS(ind_noi[s],i);
         int count1 = 0;
         for(int j = 0; j < p-1; j ++){
           if (gamma_M(ind_noi[j],i) == 1){
             Sig11(count,count1) = Sig(ind_noi[s],ind_noi[j]);
             count1 +=1;
           }
           
         }
         count +=1;
       }     
     }
     //Rcout<<"Sig11 = "<<Sig11<<std::endl;
    //Rcout<<"Sig12 = "<<Sig12<<std::endl;
     invC11 = Sig11 - Sig12 * Sig12.t()/Sig(i,i);
     
     Ci = (xS(i,i)+xlambda)*invC11+diagmat(tau_temp);
    //Rcout << Ci << std::endl;
     Ci_chol = chol(Ci);
     NumericVector rnd = rnorm(count);
     arma::mat rdn(rnd);

     arma::vec beta = -solve(Ci,S12) + solve(Ci_chol,rdn);
     count = 0;
     for(int s = 0; s < p-1; s ++){
       if(gamma_M(ind_noi[s],i) == 1){
         xOmega(ind_noi[s],i) = beta[count];
         xOmega(i, ind_noi[s]) = beta[count];
         count +=1;
       }else{
         xOmega(ind_noi[s],i) = 0;
         xOmega(i,ind_noi[s]) = 0;
       }
     }
     }
     
     
     // update Covariance matrix according to one-column change of precision matrix
     
     Sig = inv(xOmega);
    // Rcout << Sig<<std::endl;
   }
   
   //update lambda 
   b_post = b_lambda + accu(abs(xOmega))/2;
   NumericVector lam = rgamma(1,a_post, 1.0/b_post);
         
   //return List::create(Named("Omega") = xOmega, Named("Sigma") = Sig, Named("lambda") = lambda);
   return wrap(gam[0]);
                     
}