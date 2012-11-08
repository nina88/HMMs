#include "f1.h"
#include "f2.h"
#include "getcolumn.h"

SEXP f2(SEXP C){
  using namespace Rcpp;
  int c = as<int>(C);
  
  NumericVector x = f1(Rcpp::wrap(c));
 x[0] = x[0] + 1;

  
  return x;
}