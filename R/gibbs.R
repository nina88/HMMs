####################################################
#Public functions
####################################################
#' Gibbs sampling algorithm
#'
#' @useDynLib HMMs
#' @aliases class_check checkpoint_files log_prior log_likelihood initialise initialise_transition_matrices
#' @param y An hmm_fasta object
#' @param iter number of iterations
#' @param prior "prior_parameters" class object
#' @param burnin amount of burnin required
#' @param thin amount of thinning wanted
#' @param checkpoint "hmm_checkpoint" object if checkpointing wanted or NULL if not wanted
#' @keywords character
#' @export

#initialise_output

gibbs <- function(y, iter, prior, burnin, thin, checkpoint = NULL)
{
  hour=class_check(y, iter, prior, burnin, thin, checkpoint)
  
  ##### sort out joins
  if (length(y$join)==2) {
    joins=0
    } else { 
      joins=y$join[2:(length(y$join)-1)]
      } 
  
  #### find sequence and f and length of y (n)
  f=y$level
  y=y$fasta_seq
  y=factor(y, levels=1:f)
  n = length(y)
  r=prior$r
  
  ##### check for lambda and P existing/ initialise them
  init = initialise(prior, checkpoint)
  lambda = init$lambda
  P = init$P
  count=init$count
  segment1=init$segment1
  RT=init$RT
  print(RT)
  
  ### Prior parameter for lambda
  b=prior$b
  P.mat=prior$P.mat
    
  #### check not finished
  if (count>=iter+1) {
    message("Maximum iterations reached")
    break
  }
  
  #### create storage vectors and matrices
  ht = hour/thin
  posterior.temp = numeric(ht)
  lambda.store = matrix(0,nrow = ht, ncol=r^2)
  P.store = matrix(0, nrow = ht, ncol=r*f^2)
  segment.store=matrix(0, nrow = ht, ncol=n)
  store_count = 1
  
  ##### FB and dirichlets
  for (i in count:iter) {
    print(i)
    
    ##### Forward backward
    st = FB(y, lambda, P)
    segment2 = st$s
        
    #### label switching check
    if (i==1){
      segment2 = segment2
      P = P
      lambda = lambda
      RT = initialise_RT(r, segment2, n)
     
      }
    else {
      segm = label_switch(RT, segment2, n, lambda, P)
      segment2 = segm$s
      
      P = segm$P
      lambda = segm$lambda
      RT=segm$RT
     
      }
    segment2 = factor(segment2, levels = 1:r)
    y = factor(y, levels = 1:f)
    
        
    ### find parameters for dirichlets
    s.trans = table(segment2[1:(n-1)], segment2[2:n]) 
    y.trans = table(y[1:(n-1)], y[2:n], segment2[2:n]) 
        
    # take off transitions between proteins if real data
    if(all(joins!=0)){
      for (m in 1:length(joins)){
        s.trans[segment2[joins[m]], segment2[joins[m]+1]] = s.trans[segment2[joins[m]], segment2[joins[m]+1]]-1
        y.trans[y[joins[m]], y[joins[m]+1],segment2[joins[m]+1]] = y.trans[y[joins[m]], y[joins[m]+1],segment2[joins[m]+1]]-1
        }
      }
    #### simulate from dirichlets
    ### find p
    for (k in 1:r){
      for (j in 1:f){
        P[j,,k] = rdiric(1, P.mat[j,,k]+y.trans[j,,k])
        }
      }
    ### find lambda
    for (k in 1:r){
      lambda[k,] = rdiric(1, b[k,]+s.trans[k,])
      }
    
   
    if (i > burnin & i%%thin==0){
      
      ### log posterior
      log.prior = log_prior(P, lambda, P.mat, b)
      loglike.store = log_likelihood(P, lambda, y.trans, s.trans)  
      post = loglike.store+log.prior
      
      ### temp store
      posterior.temp[store_count]=post
      P.store[store_count,]=as.vector(P)
      lambda.store[store_count,]=as.vector(lambda)
      
      
      segment.store[store_count,]=segment2
        
      ## write lambda and P to file for checkpointing
      if (i%%hour==0) {
        checkpoint_files(lambda, P, segment1, posterior.temp, P.store, lambda.store, segment.store, i, r, checkpoint, RT)
      }
      store_count = store_count + 1
      if (store_count == ht+1){
        posterior.temp = numeric(ht)
        lambda.store = matrix(0,nrow = ht, ncol=r^2)
        P.store = matrix(0,nrow = ht, ncol=r*f^2)
        segment.store=matrix(0, nrow = ht, ncol=n)
        store_count=1
      }
      } 
  }
}

##Pass output files as arguments
##Restart argument =TRUE or FALSE
##Change to initialise_state

#' Initialising the prior
#'
#' @param P.mat prior value for the observed sequence
#' @param mu prior mean of the hidden sequence
#' @param s prior standard deviation for the hidden sequence
#' @return \item{prior }{"prior_parameters" object} 
#' @keywords character
#' @export
#' 
initialise_prior <- function(P.mat, mu, s)
{
  r=dim(P.mat)[3]
  f=dim(P.mat)[1]
  c = ((mu^2*(1-mu))/(s^2))-mu
  d = (c*(1-mu))/((r-1)*mu)
  b = matrix(d, ncol=r, nrow=r)
  diag(b) = c
  prior = list(P.mat=P.mat, b=b, r=r, f=f)
  class(prior) = "prior_parameters"
  message(paste("You have specified r =",r,"and f =",f))
  return(prior)
}

#' Initialising checkpointing
#'
#' @param filename name of file for checkpointing
#' @param hour how often checkpointing is required in number of iterations
#' @return \item{cp }{"hmm_checkpoint" object} 
#' @keywords character
#' @export
#'

#########
initialise_checkpoint = function(filename, hour) 
{ 
  cp=list(filename=filename, hour=hour)
  class(cp)="hmm_checkpoint"
  return(cp)
}

####################################################
#Private functions
####################################################

class_check <- function(y, iter, prior, burnin, thin, checkpoint)
{
  if (class(y)!="hmm_fasta"){
    stop("Object y not from correct class")
  }
  
  if (class(prior)!="prior_parameters"){
    stop("Object prior not from correct class")
  }
  
  if (is.null(checkpoint)){
    message("Note that you are not checkpointing")
    hour=iter
  } else if (class(checkpoint)!="hmm_checkpoint")
  {
    stop("Object checkpoint not from correct class")
  }
  
  ##### checkpoint arguments
  if (!is.null(checkpoint)){
    hour = checkpoint$hour
    expect_that(checkpoint$filename, matches(".Rdata"))
    if (hour%%thin!=0){
      stop("Hour must be a multiple of thin")
    }
    if (iter%%hour!=0){
      stop("Iter must be a multiple of hour")
    }
    if (burnin%%hour!=0){
      stop("Hour must be a multiple of burnin")
    }
    
  }
  return(hour)
}

##Store thin, burnin in checkpointing

checkpoint_files <- function(lambda, P, segment1, posterior.temp, P.store, lambda.store, segment.store, i, r, checkpoint, RT)
{
  count=i+1
  if (!is.null(checkpoint)){
    save(lambda,P,segment1,count,RT,file=checkpoint$filename)
  } 
  write.table(posterior.temp,file=paste("output",r,".csv",sep=""),append=T,row.names=F,col.names=F)
  write.table(P.store,file=paste("P.store",r,".csv",sep=""),append=T,row.names=F,col.names=F)
  write.table(lambda.store,file=paste("lambda.store",r,".csv",sep=""),append=T,row.names=F,col.names=F)
  write.table(segment.store,file=paste("segment.store",r,".csv",sep=""),append=T,row.names=F,col.names=F)  
}
##### find log prior

log_prior=function(P, lambda, P.mat, b)
{
  r = dim(lambda)[1]
  P.prior = sum((P.mat-1)*log(P)) + sum(lgamma(colSums(P.mat))-colSums(lgamma(P.mat)))
  lambda.prior = sum((b-1)*log(lambda))+r*lgamma(sum(b[1,]))-r*sum(lgamma(b[1,]))
  log.prior = P.prior+lambda.prior
  return(log.prior)
}

### find log likelihood
log_likelihood=function(P, lambda, y.trans, s.trans)
{
  r = dim(lambda)[1]
  f = dim(P)[1]
  ###### calculate log likelihood for transition probabilities (P)
  P.loglike = numeric(r)
  for (l in 1:r){
    P.loglike[l]=sum(y.trans[,,l]*log(P[,,l]))
  }
  ###### calculate log likelihood for transition probabilities (lambda)
  lambda.loglike = sum(s.trans*log(lambda))
  
  ##### find log like and log post
  loglike.store = sum(P.loglike)+lambda.loglike+log(1/(f*r))
  return(loglike.store)
}


### find log likelihood
log_likelihood_P=function(P, y.trans)
{
  f = dim(P)[1]
  r = dim(P)[3]
  print(c(r,f))
  print(P)
  print(y.trans)
  ###### calculate log likelihood for transition probabilities (P)
  P.loglike = numeric(r)
  for (l in 1:r){
    P.loglike[l]=sum(y.trans[,,l]*log(P[,,l]))
  }
  ##### find log like 
  loglike.store = sum(P.loglike)+log(1/f)
  print(loglike.store)
  return(loglike.store)
}

###########


initialise<- function(prior, checkpoint)
{ 
  r=prior$r
  f=prior$f
  if (is.null(checkpoint)) {
    message("Making files")
    transition_matrices = initialise_transition_matrices(prior)
    lambda = transition_matrices$lambda
    P = transition_matrices$P
    count = 1
    segment1 = NA
    RT=NA
  } else if (!is.null(checkpoint) & file.exists(checkpoint$filename)){
    message("Using existing files")
    load(checkpoint$filename)
  } else {
    message("Making files")
    transition_matrices = initialise_transition_matrices(prior)
    lambda = transition_matrices$lambda
    P = transition_matrices$P
    count = 1
    segment1 = NA
    RT=NA
  }
  return(list(lambda = lambda, P = P, count = count, segment1 = segment1, RT=RT))
}

#initialise_transition_matrices
initialise_transition_matrices <- function(prior)
{ 
  # P
  r=prior$r
  f=prior$f
  P.mat=prior$P.mat
  P=array(0,c(f,f,r))
  for(j in 1:r){
    for (i in 1:f){
    P[,,j] = rdiric(1,P.mat[i,,j])
    }
  }
  
  ### lambda 
  b = prior$b
  lambda=matrix(0, nrow=r, ncol=r)
  for (k in 1:r){
    lambda[k,] = rdiric(1, b[k,])
  }
  return(list(lambda=lambda, P=P))
}

