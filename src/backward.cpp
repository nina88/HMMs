#include "backward.h"

SEXP backward(SEXP A, SEXP B, SEXP myArray, SEXP myDims){
    
  arma::mat f = Rcpp::as<arma::mat>(A);
  arma::mat lambda=Rcpp::as<arma::mat>(B);
  Rcpp::NumericVector back(myArray);
  Rcpp::IntegerVector  dimsback(myDims);

  int r = f.n_rows;
  int n = f.n_cols;
  int i,j,k;

  arma::cube Back(back.begin(), dimsback[0], dimsback[1], dimsback[2], false);

  for ( i = n-2; i >= 0; i-- ) {
 for ( k = 0; k < r; k++ ) {
   for ( j = 0; j < r; j++ ) {
  Back.slice(i)(j,k)=lambda(j,k)*f(j,i);
}
Back.slice(i)(arma::span(),k)=Back.slice(i)(arma::span(),k)/sum(Back.slice(i)(arma::span(),k));
 }  
}

 return Rcpp::wrap(Back);
}

 



