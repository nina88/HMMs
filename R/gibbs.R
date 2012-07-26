#initialise_output

gibbs <- function(y, iter, prior, r, checkpoint = NULL)
{
  if (class(y)!="hmm_fasta"){
    stop("Object y not from correct class")
    }
  
  if (class(prior)!="prior_parameters"){
    stop("Object prior not from correct class")
    }
  
  if (is.null(checkpoint)){
    message("Note that you are not checkpointing, removing a burn in, or thinning")
    checkpoint=initialise_checkpoint(0,1,1)
    }
  
  if(class(checkpoint)!="hmm_checkpoint"){
    stop("Object checkpoint not from correct class")    
    }
  
  ##### checkpoint arguments
  hour = checkpoint$hour
  thin = checkpoint$thin
  burnin = checkpoint$burnin
  
  ##### sort out joins
  if (length(y$join)==2) {
    joins=0
    } else { 
      joins=y$join[2:(length(y$join)-1)]
      } 
  
  #### find sequence and f and length of y (n)
  f=y$level
  y=y$fasta_seq
  y=factor(y,levels=1:f)
  n = length(y)
  
  ##### check for lambda and P existing/ initialise them
  init = initialise(prior, f, r)
  lambda = init$lambda
  P = init$P
  count=init$count
  segment1=init$segment1
  
  ### Prior parameter for lambda
  b=prior$b
  a=prior$a
    
  #### check not finished
  if (count==iter+1) break
  
  #### create storage vectors and matrices
  print(hour)
  print(thin)
  posterior.temp = numeric(hour/thin)
  lambda.store = matrix(0,nrow = hour/thin, ncol=r^2)
  P.store = matrix(0,nrow = hour/thin, ncol=r*f^2)
  segment.store=matrix(0, nrow = hour/thin, ncol=n)
  store_count = 1
  
  #### store initial values
  #if ((count-1)==1){
    #lambda.store[1,]=as.vector(lambda)
    #P.store[1,]=as.vector(P)
    #}
  
  ##### FB and dirichlets
  for (i in count:iter) {
    print(i)
    
    ##### Forward backward
    st = FB(y, lambda, P)
    segment2 = st$s
        
    #### label switching check
    if (i==2){
      segment2 = segment2
      P = P
      lambda = lambda
      }
    else {
      segm = label_switch(segment2, segment1, r,lambda, P, f)
      segment2 = segm$segment.store2
      P = segm$P
      lambda = segm$lambda
      }
    segment2 = factor(segment2, levels = 1:r)
    y = factor(y, levels = 1:f)
    segment1 = segment2
        
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
        P[j,,k] = rdiric(1, a+y.trans[j,,k])
        }
      }
    ### find lambda
    for (k in 1:r){
      lambda[k,] = rdiric(1, b[k,]+s.trans[k,])
      }
    
   
    if (i >= burnin & i%%thin==0){
      
      ### log prior
      P.prior = sum((a-1)*log(P))+r*f*lgamma(f*a)-r*(f^2)*lgamma(a)
      lambda.prior = sum((b-1)*log(lambda))+r*lgamma(sum(b[1,]))-r*sum(lgamma(b[1,]))
      log.prior = P.prior+lambda.prior
        
      ###### calculate log likelihood for transition probabilities (P)
      P.loglike = numeric(r)
      for (l in 1:r){
        P.loglike[l]=sum(y.trans[,,l]*log(P[,,l]))
        }
      ###### calculate log likelihood for transition probabilities (lambda)
      lambda.loglike = sum(s.trans*log(lambda))
    
      ##### find log like and log post
      loglike.store = sum(P.loglike)+lambda.loglike+log(1/(f*r))
      post = loglike.store+log.prior
      
      ### temp store
      print(store_count)
      posterior.temp[store_count]=post
      print(posterior.temp)
      P.store[store_count,]=as.vector(P)
      print(P.store)
      lambda.store[store_count,]=as.vector(lambda)
      print(lambda.store)
      segment.store[store_count,]=segment2
      print(segment.store)
        
      ## write lambda and P to file for checkpointing
      if (i%%hour==0) {
        checkpoint(lambda, P,segment1, posterior.temp, P.store, lambda.store, segment.store, i, r)
      }
      store_count = store_count + 1
      if (store_count == hour/thin+1){
        store_count=1
      }
      } 
  }
    return(list(lambda = lambda, P = P))
}

##Change .txt to csv
##Pass output files as arguments
##Restart argument =TRUE or FALSE
##Change to initialise_state


initialise<- function(prior, f, r)
{
  if (file.exists("checkpoint.Rdata")==T){
    message("Using existing files")
    load("checkpoint.Rdata")
    } else {
      message("Making files")
      transition_matrices = initialise_transition_matrices(prior, r, f)
      lambda = transition_matrices$lambda
      P = transition_matrices$P
      count = 2
      segment1 = NA
      }
  return(list(lambda = lambda, P = P, count = count, segment1 = segment1))
}

#initialise_transition_matrices
initialise_transition_matrices <- function(prior, r, f)
{ 
  # P
  a=prior$a
  P=array(0,c(f,f,r))
  for(j in 1:r){
    P[,,j] = rdiric(f,rep(a,f))
  }
  
  ### lambda 
  b = prior$b
  diag(b) = c
  lambda=matrix(0, nrow=r, ncol=r)
  for (k in 1:r){
    lambda[k,] = rdiric(1, b[k,])
  }
  return(list(lambda=lambda, P=P))
}


initialise_prior <- function(a, mu, s, r)
{
  c = ((mu^2*(1-mu))/(s^2))-mu
  d = (c*(1-mu))/((r-1)*mu)
  b = matrix(d, ncol=r, nrow=r)
  diag(b) = c
  prior = list(a=a, b=b)
  class(prior) = "prior_parameters"
  return(prior)
}

##Store thin, burnin in checkpointing

checkpoint <- function(lambda, P, segment1, posterior.temp, P.store, lambda.store, segment.store, i, r)
{
  count=i+1
  save(lambda,P,segment1,count,file="checkpoint.Rdata")
  write.table(posterior.temp,file=paste("output",r,sep=""),append=T,row.names=F,col.names=F)
  write.table(P.store,file=paste("P.store",r,sep=""),append=T,row.names=F,col.names=F)
  write.table(lambda.store,file=paste("lambda.store",r,sep=""),append=T,row.names=F,col.names=F)
  write.table(segment.store,file=paste("segment.store",r,sep=""),append=T,row.names=F,col.names=F)	
}


#########
initialise_checkpoint = function(burnin, thin, hour) 
{
  if (hour%%thin!=0){
    stop("Hour must be a multiple of thin")
  }
  
  checkpoint=list(burnin=burnin, thin=thin, hour=hour)
  class(checkpoint)="hmm_checkpoint"
  return(checkpoint)
}
