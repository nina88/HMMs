#include "forward.h"

SEXP forward(SEXP A, SEXP B, SEXP C, SEXP D, SEXP myArray, SEXP myDims){
    using namespace Rcpp;
     Rcpp::NumericMatrix f(A);
     Rcpp::NumericMatrix lambda(B);
     Rcpp::NumericVector y(C);
     Rcpp::NumericVector g(D);
     Rcpp::NumericVector P(myArray);
     Rcpp::IntegerVector  dimsP(myDims);
    
     int r = f.nrow();
     int n = f.ncol();
     int i,j,k,l;
     double a;
     Rcpp::NumericVector colvec1(r);
     Rcpp::NumericVector colvec2(r);
     Rcpp::NumericVector colvec(r);
     Rcpp::NumericVector colvec3(r);
     double s=0;
     double s2=0;
     double lml2=0;

     arma::cube p(P.begin(), dimsP[0], dimsP[1], dimsP[2], false);

  for ( i = 1; i < n; i++ ) {
   for ( k = 0; k < r; k++ ) {
       for ( j = 0; j < r; j++ ) {
   
          colvec1(j)=lambda(j,k);
          colvec2(j)=f(j,i-1);
          colvec(j)=colvec1(j)*colvec2(j);
          s=s+colvec(j);
        }
       f(k,i)=p(y(i-1)-1,y(i)-1,k)*s;
       s=0;
       
   } 

for (l=0; l<r; l++){
  colvec3(l)=f(l,i);
  s2=s2+colvec3(l);
  
}

f(_,i)=f(_,i)/s2;
lml2=lml2+log(s2);
s2=0;
   }
 
    
     return Rcpp::List::create(f,lml2); 
}


    
   
     

