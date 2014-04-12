// sample Omega columnwise, together with lambda

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP sample_Omega_W(NumericMatrix xS, IntegerMatrix ind_noi_all, double xlambda, int xn, IntegerMatrix gamma_M, NumericMatrix tau, NumericMatrix xSig, 
                  NumericMatrix Omega, double a_post, double b_lambda) {
   int p = xS.ncol();
   IntegerVector ind_noi(p-1);
   arma::mat Sig(REAL(xSig), xSig.nrow(), xSig.ncol(), true, true);
   arma::mat xOmega(REAL(Omega), Omega.nrow(), Omega.ncol(), true, true);
    
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
     //Rcout << "count = "<<count;
     if(count > 0){    
     arma::mat Sig11(p-1, p-1);
     arma::mat Sig12(p-1,1);
     arma::mat invC11(p-1,p-1);
     arma::mat Ci(p-1, p-1);
     arma::mat Ci_chol(count, count);
     arma::colvec tau_temp(p-1,1);
     arma::mat S12(count,1);
     

     for(int s = 0; s < p-1; s++){
         tau_temp[s] = tau(ind_noi[s],i);
         Sig12[s] = Sig(ind_noi[s],i);
         
         for(int j = 0; j < p-1; j ++){
             Sig11(s,j) = Sig(ind_noi[s],ind_noi[j]);
         }
     }
     Rcout << "i = "<<i <<std::endl;
     Rcout << "Sigma = "<<Sig<<std::endl;
     Rcout<<"Sig11 = "<<Sig11<<std::endl;
    Rcout<<"Sig12 = "<<Sig12<<std::endl;
    Rcout << "tau_temp = " << tau_temp <<std::endl;
    
    invC11 = Sig11 - Sig12 * Sig12.t()/Sig(i,i);
    
     Rcout << "invC11 = " << invC11 <<std::endl;
     Ci = (xS(i,i)+xlambda)*invC11+diagmat(1./tau_temp);
     Rcout << Ci << std::endl;
     
     arma::mat Ci_sub(count, count);

    count = 0;
     for(int s = 0; s < p-1; s++){
       if(gamma_M(ind_noi[s],i) == 1){
         S12[count] = xS(ind_noi[s],i);
         int count1 = 0;
         for(int j = 0; j < p-1; j ++){
           if (gamma_M(ind_noi[j],i) == 1){
             Ci_sub(count,count1) = Ci(s,j);
             count1 +=1;
           }
         }
         
         count +=1;
       }     
     }
     Rcout << "Ci_sub = "<<Ci_sub <<std::endl;
     
     Ci_chol = chol(Ci_sub);
     Rcout << "Ci_chol = "<<Ci_chol <<std::endl;
     NumericVector rnd = rnorm(count);
     arma::mat rdn(rnd);
     

     arma::vec beta = -solve(Ci_sub,S12) + solve(Ci_chol,rdn);
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
     }else{
       for(int s = 0; s < p-1; s ++){
       
         xOmega(ind_noi[s],i) = 0;
         xOmega(i,ind_noi[s]) = 0;
     
       }
     }
     
    Rcout <<"Omega = "<< xOmega <<std::endl;
    
    // update Covariance matrix according to one-column change of precision matrix
     Sig = inv(xOmega);
   }
   
   //update lambda 
   b_post = b_lambda + accu(abs(xOmega))/2;
   NumericVector lam = rgamma(1,a_post, 1.0/b_post);
         
    return List::create(Named("Omega") = xOmega, Named("Sigma") = Sig, Named("lambda") = lam[0]);
                     
}