equil <-
function(P)
{
  ev = eigen(t(P))$vectors[,1]
  return(Re(ev/sum(ev)))
}

rdiric <-
function(n,a) {
        p <- length(a)
        m <- matrix(nrow=n,ncol=p)
        for (i in 1:p) {
                m[,i] <- rgamma(n,a[i])
        }
        sumvec <- m %*% rep(1,p)
        m / as.vector(sumvec)
}

hmm.sim <-
function(n,lambda,P,f)
{
  
  ## first simulate the hidden states (segmentation)
  r = dim(lambda)[1]
  pi.lam = equil(lambda)
  s = numeric(n)
  s[1] = sample(1:r,1,replace=TRUE,prob=pi.lam)
  for(t in 2:n){
    s[t] = sample(1:r,1,replace=TRUE,prob=lambda[match(s[t-1],1:r),])
  }
  ## then simulate the observed states (conditional on s)
  pi.P = array(0,c(f,r))
  for(i in 1:r){
    pi.P[,i] = equil(P[,,i])
  }
  y = numeric(n)
  y[1] <- sample(1:f,1,replace=TRUE,prob=pi.P[,s[1]])
  for(t in 2:n){
    y[t] <- sample(1:f,1,replace=TRUE,prob=P[y[t-1],,s[t]])
  }
  list(s=s,y=y)
}
