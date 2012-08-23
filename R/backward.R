require(inline)
require(RcppArmadillo)

src <- '

Rcpp::NumericMatrix f(A);
Rcpp::NumericMatrix lambda(B);
Rcpp::NumericVector back(myArray);
Rcpp::IntegerVector  dimsback(myDims);

int r = f.nrow();
int n = f.ncol();
int i,j,k;
NumericMatrix vec(r,n);

arma::cube Back(back.begin(), dimsback[0], dimsback[1], dimsback[2], false);

 for ( i = n-2; i >= 0; i-- ) {
    for ( k = 0; k < r; k++ ) {
      for (j=0; j<r; j++){
         Back(j,k,i)=lambda(j,k)*f(j,i);  
         vec(k,i)+=Back(j,k,i);   
         } 
     }
  }


 for ( i = n-2; i >= 0; i-- ) {
    for ( k = 0; k < r; k++ ) {
      for (j=0; j<r; j++){
         Back(j,k,i)/=vec(k,i);    
         } 
     }
  }

return Rcpp::wrap(Back);
'


backward <- cxxfunction(signature(A = "numeric",B="numeric",myArray="numeric",myDims="integer"), body = src, plugin="RcppArmadillo")
#backward(f,lambda,back,dim(back))
