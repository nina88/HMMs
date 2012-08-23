####################################################
#Public functions
####################################################
#' The forward--backward algorithm for the power posterior method
#'
#' @param y An hmm_fasta object
#' @param lambda hidden sequence transition matrix
#' @param P array of transition matrices for observed sequence
#' @param t temperature parameter
#' @return \item{s }{segmentation} 
#'  \item{xi }{normalising constant}
#'  \item{back}{backwards probabilities}
#' @keywords character
#' @export
##### took out lambdas to the power of t
FBpower = function(y, lambda, P, t)
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
 
  #####################################################################

  ## initialise the forward probabilities 
  pi.P = array(0,c(h,r))
  
  for(i in 1:r){
      pi.P[,i] = equil(P[,,i])
      }
  
  pi.P=pi.P^t
  P=P^t
  pi.lambda = equil(lambda)
  
  xi[1] = sum(pi.P[y[1], ]*pi.lambda)
  f[ ,1] = (pi.P[y[1], ]*pi.lambda)/xi[1] 
  
  #####################################################################
  
  ## calculate the forward probabilities using forward recursions

  for(i in 2:n){
    for(k in 1:r){
      f[k, i] = P[y[i-1], y[i], k]*sum(lambda[ ,k]*f[ ,i-1])
    } 
    xi[i]=sum(f[ ,i])
    f[ ,i] = f[ ,i]/xi[i]
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
    s[i] = sample(1:r, 1, prob = back[,s[i+1],i])
  }
  ## return the global estimate
  return(list(s = s, xi=xi, back=back))
}

####################################
#' The power posterior method
#'
#' @aliases checkpoint_files_power.R class_check_power.R initialise_power.R initialise_transition_matrices_power.R
#' @param N number of temperature parameters required minus 1
#' @param y An hmm_fasta object
#' @param prior "prior_parameters" class object
#' @param m number of iterations per temperature parameter
#' @param r number of segment types
#' @param burnin amount of burn in required
#' @param checkpoint "hmm_checkpoint" object or NULL if checkpointing not required
#' @return \item{logpPP }{marginal likelihood for r} 
#' @keywords character
#' @export
#' 
powerpost=function(N, prior, m, r, y, burnin, checkpoint=NULL){
  
  ### insert class checks here  
  class_check_power(y, prior, checkpoint)
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
  init = initialise_power(prior, f, r, checkpoint)
  lambda = init$lambda
  P = init$P
  expectation=init$expectation
  count=init$count
  
  
  ### Prior parameters
  b=prior$b
  a=prior$a
     
 ### step 2: a: set up temperature parameter t_i (in loop)
  loglike.store=numeric(m)
  
 ### step 2: b: generate sample of theta from the power posterior
  for (i in count:N){
        t=(i/N)^4
        print(t)
        #### need a,b s.trans and y.trans from gibbs sampling
        ### storage for P and lambda
        lambda.store=matrix(0,nrow=m, ncol=r^2)
        P.store=matrix(0,nrow=m, ncol=r*f^2)
        
        for (h in 1:m){  
            st = FBpower(y, lambda, P,t)
            segment2 = st$s
            segment2 = factor(segment2, levels = 1:r)
            y = factor(y, levels = 1:f)
            ### find parameters for dirichlets
            s.trans = table(segment2[1:(n-1)], segment2[2:n]) 
            y.trans = table(y[1:(n-1)], y[2:n], segment2[2:n]) 
   
            # take off transitions between proteins
            for (q in 1:length(joins)){
              s.trans[segment2[joins[q]], segment2[joins[q]+1]] = s.trans[segment2[joins[q]], segment2[joins[q]+1]]-1
              y.trans[y[joins[q]], y[joins[q]+1],segment2[joins[q]+1]] = y.trans[y[joins[q]], y[joins[q]+1],segment2[joins[q]+1]]-1
                    }
    
            ### find p
            for (k in 1:r){
                for (j in 1:f){
                    P[j,,k] = rdiric(1, a+t*y.trans[j,,k])
                }  
            }
            ### find lambda
            for (k in 1:r){
                lambda[k,] = rdiric(1, b[k,]+t*s.trans[k,])
            } 
            lambda.store[h,]=as.vector(lambda)
            P.store[h,]=as.vector(P)
        
            ### step 2: c: estimate the expectation
            ### calculate log likelihood
            loglike = log_likelihood(P, lambda, y.trans, s.trans, r, f)
            loglike.store[h]=loglike
          
        } 
        ### step 2 c (not need be in l)
        #write.table(loglike.store,file="loglike", append=T, row.names=F, col.names=F)
        expectation[i+1]=mean(loglike.store[burnin:m])
        
        ### step 2 d  new lambda and P
        expectation.lambda=apply(lambda.store[burnin:m,],2,mean)
        expectation.P=apply(P.store[burnin:m,],2,mean)
        lambda=matrix(expectation.lambda,nrow=r,ncol=r)
        P=array(expectation.P,c(f,f,r))   
        
        ##### checkpointing system
        checkpoint_files_power(lambda, P, expectation, count, i, checkpoint)
  }
  
  #### step 3 calculate logp_PP(y|r)
  T=numeric(N+1)
  for (i in 0:N){
    T[i+1]=(i/N)^4
  }
  logpPP.store=numeric(N)
  for (i in 1:N){
      logpPP.store[i]=((T[i+1]-T[i])*((expectation[i+1]+expectation[i])/2))
      }
  logpPP=sum(logpPP.store)
  write.csv(logpPP,file="logpPP.txt")
  return(list(logpPP=logpPP))
}

#########
#' Initialising checkpointing for the power posterior method
#'
#' @param filename filename
#' @return \item{cp }{"hmm_checkpoint" object} 
#' @keywords character
#' @export
#'
initialise_checkpoint_power = function(filename) 
{ 
  cp=list(filename=filename)
  class(cp)="hmm_checkpoint"
  return(cp)
}

####################################################
#Private functions
####################################################

####################################################################################################
####################################################################################################
### note needs to be changed

checkpoint_files_power <- function(lambda, P, expectation, count, i, checkpoint)
{
  count=i+1
  if (!is.null(checkpoint)){
    save(lambda, P, expectation, count, file=checkpoint$filename)
  }  
}

class_check_power <- function(y, prior, checkpoint)
{
  if (class(y)!="hmm_fasta"){
    stop("Object y not from correct class")
  }
  
  if (class(prior)!="prior_parameters"){
    stop("Object prior not from correct class")
  }
  
  if (is.null(checkpoint)){
    message("Note that you are not checkpointing")
  } else if (class(checkpoint)!="hmm_checkpoint")
  {
    stop("Object checkpoint not from correct class")
  }
  
  ##### checkpoint arguments
  if (!is.null(checkpoint)){
    expect_that(checkpoint$filename, matches(".Rdata"))
  }
}

initialise_power<- function(prior, f, r, checkpoint, N)
{ 
  if (is.null(checkpoint)) {
    message("Making files")
    transition_matrices = initialise_transition_matrices_power(r, f)
    lambda = transition_matrices$lambda
    P = transition_matrices$P
    count = 0
    expectation=numeric(N+1)

  } else if (!is.null(checkpoint) & file.exists(checkpoint$filename)){
    message("Using existing files")
    load(checkpoint$filename)
  } else {
    message("Making files")
    transition_matrices = initialise_transition_matrices_power(r, f)
    lambda = transition_matrices$lambda
    P = transition_matrices$P
    count = 0
    expectation=numeric(N+1)
  }
  return(list(lambda = lambda, P = P, count = count, expectation = expectation))
}


#initialise_transition_matrices
initialise_transition_matrices_power <- function(r, f)
{ 
  P=array(1/f,c(f,f,r)) 
  diag.prob = 0.9    ## probability of staying in state i
  lambda = matrix((1-diag.prob)/(r-1),nrow=r,ncol=r)
  for(k in 1:r){
    lambda[k,k]=diag.prob
  }
  return(list(lambda=lambda, P=P))
}

#########

