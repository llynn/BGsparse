#include <Rcpp.h> 
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector colMin_Max(NumericMatrix X, int c) {
  if ( c == 1){ // max
    int ncol = X.ncol();
    Rcpp::NumericVector out(ncol);
    for (int col = 0; col < ncol; col++){
        out[col]=Rcpp::max(X(_, col)); 
    } 
    return wrap(out);
    
  }else{ // min
    int ncol = X.ncol();
    Rcpp::NumericVector out(ncol);
    for (int col = 0; col < ncol; col++){
        out[col]=Rcpp::min(X(_, col)); 
    } 
    return wrap(out);
  } 
  
} 