// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double determinant(NumericMatrix x) {
  arma::mat m(REAL(x), x.nrow(), x.ncol(), false, true);
  double output = det(m);
  return output;
}

Function R_det("det");

// [[Rcpp::export]]
double R_det_(NumericMatrix x) {
  return as<double>(R_det(x));
}