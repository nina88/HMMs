####################################################
#Public functions
####################################################
#' The forward--backward algorithm
#'
#' @param y An hmm_fasta object
#' @param lambda hidden sequence transition matrix
#' @param P array of transition matrices for observed sequence
#' @return \item{s }{segmentation} 
#'  \item{xi }{normalising constant}
#'  \item{back}{backwards probabilities}
#' @keywords character
#' @export


FB <-
function(y, lambda, P)
{
  n = length(y)
  r = dim(lambda)[1]
  h = dim(P)[1]
  
  ## filtered probabilities
  f = matrix(0, nrow=r, ncol=n)
  
  ##### segmentation
  s = numeric(n)
  
   ## array of backward probabilities
  back = array(0 ,c(r, r, n))
  
  #### normalising constant
  xi = numeric(n)  

  # observed data loglikelihood
  lml = 0
 
  #####################################################################

  ## initialise the forward probabilities 
  
  pi.P = array(0,c(h,r))
  for(i in 1:r){
      pi.P[,i] = equil(P[,,i])
      }
  pi.lambda = equil(lambda)
  xi[1] = sum(pi.P[y[1], ]*pi.lambda)
  f[ ,1] = (pi.P[y[1], ]*pi.lambda)/xi[1]
  lml = lml + log(xi[1])
  
  #####################################################################
  
  ## calculate the forward probabilities using forward recursions

  for(i in 2:n){
    for(k in 1:r){
      f[k, i] = P[y[i-1], y[i], k]*sum(lambda[ ,k]*f[ ,i-1])
    } 
    xi[i]=sum(f[ ,i])
    f[ ,i] = f[ ,i]/xi[i]
    lml = lml +log(xi[i])
  }
  

  #####################################################################
  ## simulate backwards probabilities
  for(i in (n-1):1){
    back[, ,i] = lambda*f[,i]
    for(k in 1:r){
      back[ ,k,i] = back[ ,k,i]/sum(back[ ,k,i])
      }
  }
  
  #####################################################################
  ## Steps 5 and 6
  ## simulate segmentation 
  
  s[n] = sample(1:r, 1, prob = f[,n])
  for (i in (n-1):1){
    s[i] = sample(1:r ,1, prob = back[,s[i+1],i])
  }
  
  ## return the global estimate
  return(list(s = s, xi=xi, back=back, f=f))
}
