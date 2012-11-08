#include "f1.h"

SEXP f1(SEXP D){
  using namespace Rcpp ;
  int cons = as<int>(D);
  NumericVector x(1);
  x = cons;
  return x;
}