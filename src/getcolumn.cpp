#include "getcolumn.h"

SEXP getcolumn(SEXP m, SEXP n){
    using namespace Rcpp;
     Rcpp::NumericMatrix mat(m);
     int index = as<int>(n) ; 
     int p = mat.nrow();
     Rcpp::NumericVector colvec(p);
     int i;

     for ( i = 0; i < p; i++ ) {
     colvec(i)=mat(i,index);
     }
  
    
     return Rcpp::wrap(colvec); 
}


    
   
     
