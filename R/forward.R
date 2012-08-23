#require(inline)
#require(RcppArmadillo)

  src <- '
     Rcpp::NumericMatrix f(A);
     Rcpp::NumericMatrix lambda(B);
     Rcpp::NumericVector y(C);
     Rcpp::NumericVector g(D);
     Rcpp::NumericVector P(myArray);
     Rcpp::IntegerVector  dimsP(myDims);
    
     int r = f.nrow();
     int n = f.ncol();
     int i,j,k;

     arma::cube p(P.begin(), dimsP[0], dimsP[1], dimsP[2], false);

  for ( i = 1; i < n; i++ ) {
   for ( k = 0; k < r; k++ ) {
     for (j=0; j<r; j++){
       f(k,i)=f(k,i)+lambda(j,k)*f(j,i-1);    
     } 
     
     f(k,i)*=p(y(i-1)-1,y(i)-1,k); 
     g(i)+=f(k,i); 
   }
 }
 
    for ( i = 1; i < n; i++ ) {
    for ( k = 0; k < r; k++ ) {
     
      f(k,i)/=g(i); 
      }
    }
 
     return Rcpp::List::create(f,g);
 '

 
forward <- cxxfunction(signature(A = "numeric",B="numeric",C="numeric",D="numeric",myArray="numeric",myDims="integer"), body = src, plugin="RcppArmadillo")
#forward(f, lambda, y, xi, P, dim(P))
