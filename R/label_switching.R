####################################################
#Public functions
####################################################
#' Label switching
#'
#' @aliases initialise_RT permutations permut
#' @param RT counts for each simulated sequnce
#' @param s current segmentation sequence
#' @param r number of segment types
#' @param n length of sequence
#' @param lambda transition matrix for hidden sequence
#' @param P array of transition matrices for observed sequence
#' @param f number of possible states of observed sequence
#' @return \item{RT}{counts}
#' \item{s}{segmentation after label switching corrected for} 
#'\item{P }{tranition matrix P after label switching corrected for }
#'  \item{lambda}{trabsition matric lambda after label switching corrected for}
#' @keywords character
#' @export



label_switch=function(RT, s, r, n, lambda, P, f)
{
  ### calculate all permutations of s and calculate sum
  vec=numeric(factorial(r))
  for (k in 1:factorial(r)){
    ##### calculate all permutations
    perms=permut(r, s, k)
    for (i in 1:n){
      vec[k]=vec[k]+RT[perms[i],i]
    }  
  }
  
  #### update RT
  d=which.max(vec)
  perms=permut(r, s, d)
  
  for (i in 1:n){
    for (j in 1:r){
      if (perms[i]==j){
        RT[j,i]=RT[j,i]+1
      }
    }
  }
  
  z1=order_of_firsts(s,r)
  z2=order_of_firsts(perms,r)
  new.P=array(0, c(f,f,r))
  new.lambda=matrix(0, nrow=r, ncol=r)
  if (identical(s,perms)){
    new.P=P
    new.lambda=lambda
  } else {
    for (k in 1:r){
      new.P[,,z2[k]]=P[,,z1[k]]
      for (l in 1:r){
        new.lambda[z2[k],z2[l]]=lambda[z1[k],z1[l]]
      }
    }
  }
  
  return(list(RT=RT,s=perms,lambda=new.lambda,P=new.P))
}



####################################################
#Private functions
####################################################
initialise_RT=function(r,s,n)
{
  RT=matrix(0,nrow=r,ncol=n)
  
  ### update RT for inital sequence
  for (i in 1:n){
    for (j in 1:r){
      if (s[i]==j){
        RT[j,i]=RT[j,i]+1
      }
    }
  }
  return(RT)
}

permut <-
  function(r, segments, i) {
    
    perm = permutations(r,r)
    vec=numeric(length(segments))
    new_vec = segments
    
    for(j in 1:r) {
      new_vec[segments == j] = perm[i,j]
    }
    return(new_vec) 
  }


##### this function is taken from the package gregmisc to avoid dependencies

permutations <-
  function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
  {
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
          0) 
      stop("bad value of n")
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
          0) 
      stop("bad value of r")
    if (!is.atomic(v) || length(v) < n) 
      stop("v is either non-atomic or too short")
    if ((r > n) & repeats.allowed == FALSE) 
      stop("r > n and repeats.allowed=FALSE")
    if (set) {
      v <- unique(sort(v))
      if (length(v) < n) 
        stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) 
      sub <- function(n, r, v) {
        if (r == 1) 
          matrix(v, n, 1)
        else if (n == 1) 
          matrix(v, 1, r)
        else {
          inner <- Recall(n, r - 1, v)
          cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                    ncol = ncol(inner), nrow = nrow(inner) * n, 
                                                    byrow = TRUE))
        }
      }
    else sub <- function(n, r, v) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        X <- NULL
        for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                          1, r - 1, v[-i])))
        X
      }
    }
    sub(n, r, v[1:n])
  }



