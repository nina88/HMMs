####################################################
#Public functions
####################################################
#' Label switching
#'
#' @aliases order_of_firsts.R permutations.R permut.R
#' @param segment.store2 segmentation at iteration i
#' @param segment.store1 segmentation at iteration i-1
#' @param r number of segment types
#' @param lambda transition matrix for hidden sequence
#' @param P array of transition matrices for observed sequence
#' @param f number of possible states of observed sequence
#' @return \item{segment.store2 }{segmentation after label switching corrected for} 
#'  \item{P }{tranition matrix P after label switching corrected for }
#'  \item{lambda}{trabsition matric lambda after label switching corrected for}
#' @keywords character
#' @export

label_switch <-
  function(segment.store2,segment.store1, r, lambda, P, f)
  {
    
    ##### vec for matches, matrices for storage
    x=numeric(factorial(r))
    new.P=array(0, c(f,f,r))
    new.lambda=matrix(0, nrow=r, ncol=r)
    
    ### segment.store2 is i and 1 is i-1
    for (i in 1:factorial(r)){
      ##### calculate all permutations
      perms=permut(r, segment.store2, i)
      x[i]=sum(perms == segment.store1)
    }
    
    ### find location of maximum in x
    y=which.max(x)
    
    #### return correct labelling
    segment.store3=permut(r, segment.store2, y)
    
    z1=order_of_firsts(segment.store2,r)
    z2=order_of_firsts(segment.store3,r)
    
    if (identical(segment.store2,segment.store3)){
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
    
    return(list(segment.store2=segment.store3, P=new.P, lambda=new.lambda))
    
  }

####################################################
#Private functions
####################################################

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

##### 

order_of_firsts <-
  function(data, length){
    first = array()
    for(i in 1:length){
      first[i] = which(data==i)[1]
    }
    order(first)
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
