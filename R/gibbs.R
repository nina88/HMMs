#initialise_prior
#initialise_output
#initialise_checkpoint = function(hour) {
    #do something
  #  class("hmm_checkpoint")
  #  return(somethiong)
#}

#ic  = initialise_checkpoint(10)
#gibbs(y, iter, ic)
gibbs <- function(y, iter, hour, prior, r)##, checkpoint =NULL, prior)
{
    if (class(y)!="hmm_fasta"){
        stop("Object y not from correct class")
    }
    
    if (class(prior)!="prior_parameters"){
      stop("Object prior not from correct class")
    }
    
    if(!is.null(checkpoint)){
        
    }
    f=y$level
    if (length(y$join)==2) {
        joins=0
    } else { 
        joins=y$join[2:(length(y$join)-1)]
    } 
    y=y$fasta_seq
    y=factor(y,levels=1:f)
    
    ##### check for lambda and P existing/ initialise them
    init = initialise(prior, f, r)
    lambda = init$lambda
    P = init$P
    count=init$count
    segment1=init$segment1
    ### Prior parameter for lambda
    b=init$b
    a=prior$a
    #### check not finished
    if (count==iter+1) break
    n = length(y)
    posterior.temp=numeric(hour)
    lambda.store=matrix(0,nrow=hour, ncol=r^2)
    P.store=matrix(0,nrow=hour, ncol=r*f^2)
    segment.store=matrix(0, nrow=hour, ncol=n)
    if ((count-1)==1){
        lambda.store[1,]=as.vector(lambda)
        P.store[1,]=as.vector(P)
    }
    
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
        if (i%%hour==0){
            posterior.temp[hour]=post
            P.store[hour,]=as.vector(P)
            lambda.store[hour,]=as.vector(lambda)
            segment.store[hour,]=segment2
        } else{
            posterior.temp[i%%hour]=post
            P.store[i%%hour,]=as.vector(P)
            lambda.store[i%%hour,]=as.vector(lambda)
            segment.store[i%%hour,]=segment2
        }
        
        ## write lambda and P to file for checkpointing
        if (i%%hour==0){
            checkpoint(lambda, P,segment1, posterior.temp, P.store, lambda.store, segment.store, i, r)
        }
    }  
    return(list(lambda = lambda, P = P))
}

##Change .txt to csv
##Pass output files as arguments
##Restart argument =TRUE or FALSE
##Pass prior arguments in.
##Change to initialise_state
##Have initialise_checkpoint

initialise<- function(prior, f, r)
{
    if (file.exists("lambda.txt")==T & 
        file.exists("P.txt")==T & 
        file.exists("count.txt")==T & 
        file.exists("segment1.txt")==T){
            #message(Using existing files)
            lambda = matrix(read.csv(file="lambda.txt")[,2],nrow=r,ncol=r)
            P = array(read.csv("P.txt")[,2],c(f,f,r))
            count = read.csv("count.txt")[,2]
            segment1 = read.csv("segment1.txt")[,2]
    } else {
        #message(Making files)
        transition_matrices = initialise_transition_matrices(prior, f, r)
        lambda = transition_matrices$lambda
        P = transition_matrices$P
        b = transition_matrices$b
        count = 2
        segment1 = NA
    }
    return(list(lambda = lambda, P = P, count = count, segment1 = segment1, b=b))
}

#initialise_transition_matrices
initialise_transition_matrices <- function(prior, r, f)
{ a=prior$a
  P=array(0,c(f,f,r))
  for(j in 1:r){
    P[,,j] = rdiric(f,rep(a,f))
  }
  ### Prior parameters
  mu=prior$mu
  s=prior$s
  c = ((mu^2*(1-mu))/(s^2))-mu
  d = (c*(1-mu))/((r-1)*mu)
  b = matrix(d, ncol=r, nrow=r)
  diag(b) = c
  lambda=matrix(0, nrow=r, ncol=r)
  for (k in 1:r){
    lambda[k,] = rdiric(1, b[k,])
  }
  return(list(lambda=lambda, P=P, b=b))
}


initialise_prior <- function(a, mu, s)
{
  prior = list(a=a, mu=mu, s=s)
  class(prior) = "prior_parameters"
  return(prior)
}



##Create Checkpoint.RData file using load/save
##Include everything you need to checkpoint
##Store thin, burnin

checkpoint <- function(lambda, P, segment1, posterior.temp, P.store, lambda.store, segment.store, i, r)
{
    write.csv(as.vector(lambda),file="lambda.txt")
    write.csv(as.vector(P),file="P.txt")
    write.csv(segment1,file="segment1.txt")
    write.table(posterior.temp,file=paste("output",r,sep=""),append=T,row.names=F,col.names=F)
    write.table(P.store,file=paste("P.store",r,sep=""),append=T,row.names=F,col.names=F)
    write.table(lambda.store,file=paste("lambda.store",r,sep=""),append=T,row.names=F,col.names=F)
    write.table(segment.store,file=paste("segment.store",r,sep=""),append=T,row.names=F,col.names=F)
    write.csv(i+1,file="count.txt")		
}
