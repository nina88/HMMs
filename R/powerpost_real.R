##### took out lambdas to the power of t
FBpower = function(y, lambda, p, t)
{
  n = length(y)
  r = dim(lambda)[1]
  h = dim(p)[1]
 
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
      pi.P[,i] = equil(p[,,i])
      }
  
  pi.P=pi.P^t
  p=p^t
  pi.lambda = equil(lambda)
  
  xi[1] = sum(pi.P[y[1], ]*pi.lambda)
  f[ ,1] = (pi.P[y[1], ]*pi.lambda)/xi[1] 
  
  #####################################################################
  
  ## calculate the forward probabilities using forward recursions

  for(i in 2:n){
    for(k in 1:r){
      f[k, i] = p[y[i-1], y[i], k]*sum(lambda[ ,k]*f[ ,i-1])
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

powerpost=function(N, mu, s, m, r, f, y, K, joins){
  
  expectation=numeric(N+1)
  n=length(y)
  if (file.exists("count.txt")==T & file.exists("lambda.txt")==T & file.exists("P.txt")==T & file.exists("expectation.txt")==T){
    count=read.csv("count.txt")[,2]
    lambda=matrix(read.csv(file="lambda.txt")[,2],nrow=r,ncol=r)
    P=array(read.csv("P.txt")[,2],c(f,f,r))
    expectation[1:count]=read.table(file=paste("expectation.txt"))[1:count,]
  } else {
    count=0
    #### step 1: initialise theta at prior mean
    P=array(1/f,c(f,f,r)) 
    diag.prob = 0.9    ## probability of staying in state i
    lambda = matrix((1-diag.prob)/(r-1),nrow=r,ncol=r)
    for(k in 1:r){
    lambda[k,k]=diag.prob
    }
  }
  
  ### Prior parameters
  a = 1
  c = ((mu^2*(1-mu))/(s^2))-mu
  d = (c*(1-mu))/((r-1)*mu)
  b = matrix(d, ncol=r, nrow=r)
  for (k in 1:r){
      b[k,k] = c
  }
     
 ### step 2: a: set up temperature parameter t_i (in loop)
  loglike.store=numeric(m)
  
 ### step 2: b: generate sample of theta from the power posterior
  for (i in count:N){
        #t=T[i+1]
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
            ###### calculate log likelihood for transition probabilities (P)
            P.loglike = numeric(r)
            for (l in 1:r){
                P.loglike[l]=sum(y.trans[,,l]*log(P[,,l]))
            }
  
            ###### calculate log likelihood for transition probabilities (lambda)
            lambda.loglike = sum(s.trans*log(lambda))
     
            ##### find log like and log post
            loglike = sum(P.loglike)+lambda.loglike+log(1/(f*r))
            loglike.store[h]=loglike
          
        } 
        ### step 2 c (not need be in l)
        write.table(loglike.store,file="loglike", append=T, row.names=F, col.names=F)
        expectation[i+1]=mean(loglike.store[K:m])
        
        ### step 2 d  new lambda and P
        expectation.lambda=apply(lambda.store[K:m,],2,mean)
        expectation.P=apply(P.store[K:m,],2,mean)
        lambda=matrix(expectation.lambda,nrow=r,ncol=r)
        P=array(expectation.P,c(f,f,r))   
        
        ##### checkpointing system
        write.csv(as.vector(lambda),file="lambda.txt")
        write.csv(as.vector(P),file="P.txt")
        write.table(expectation[i+1],file="expectation.txt", append=T, row.names=F, col.names=F)
        write.csv(i+1, file="count.txt")
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
  return(list(logpPP,loglike.store))
}

####################################################################################################
####################################################################################################
###other required funtions
rdiric <- function(n,a) {
        p <- length(a)
        m <- matrix(nrow=n,ncol=p)
        for (i in 1:p) {
                m[,i] <- rgamma(n,a[i])
        }
        sumvec <- m %*% rep(1,p)
        m / as.vector(sumvec)
}


equil = function(P)
{
  ev = eigen(t(P))$vectors[,1]
  return(Re(ev/sum(ev)))
}

#hmm.sim = function(n,lambda,P,f)
#{
  
  ## first simulate the hidden states (segmentation)
#  r = dim(lambda)[1]
 # pi.lam = equil(lambda)
#  s = numeric(n)
#  s[1] = sample(1:r,1,replace=TRUE,prob=pi.lam)
#  for(t in 2:n){
 #   s[t] = sample(1:r,1,replace=TRUE,prob=lambda[match(s[t-1],1:r),])
 # }
  ## then simulate the observed states (conditional on s)
#  pi.P = array(0,c(f,r))
 # for(i in 1:r){
#    pi.P[,i] = equil(P[,,i])
#  }
#  y = numeric(n)
 # y[1] <- sample(1:f,1,replace=TRUE,prob=pi.P[,s[1]])
#  for(t in 2:n){
#    y[t] <- sample(1:f,1,replace=TRUE,prob=P[y[t-1],,s[t]])
 # }
#  list(s=s,y=y)
#}
################################################################################
read.FASTA <- function(filename="")
{
  strsplit(paste(scan(filename,skip=1,what="character",comment.char=";"),collapse=""),"")[[1]]
}

## convert amino acids to hydrophilic or hydrophobic
convert4 <- function(x,frm=c(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)]),to=c(1,1,2,2,1,1,2,1,2,1,1,2,1,2,2,2,2,1,1,1))
{
  to[match(x,frm)]
}
#### convert to charge
convert5 <- function(x,frm=c(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)]),to=c(1,1,3,3,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1))
{
  to[match(x,frm)]
}
### covvery charge and hydrophobicity
convert6 <- function(x,frm=c(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)]),to=c(1,1,4,4,1,1,2,1,3,1,1,2,1,2,3,2,2,1,1,1))
{
  to[match(x,frm)]
}

##############################################################################
x=read.csv(file="data.dat")
r=x[1,2]
iterations=x[2,2]
f=x[3,2]
concat=x[4,2]

x1=read.FASTA(filename="p53.txt")
x2=read.FASTA(filename="mdm2.txt")

if (f==2){
  x1=convert4(x1)
  x2=convert4(x2)

}
if (f==3) {
x1=convert5(x1)
x2=convert5(x2)

}
if (f==4){
  x1=convert6(x1)
  x2=convert6(x2)

}

y1=c(x1,x2)
y=rep(y1,concat)
y=factor(y,levels=1:f)

###### calculate joins by arithmetic progression
#a=seq(1,(2*concat-1),2)
#b=seq(2,2*concat,2)
#joins=numeric(2*concat)
#p=numeric(concat)
#q=numeric(concat)
#p[1]=length(x1)
#q[1]=length(x2)

#if (length(p)>1){
#for (j in 2:concat){
#p[j]=p[1]+(j-1)*length(y1)
#q[j]=q[1]+(j-1)*length(y1)
#}
#}

#joins[a]=p
#joins[b]=q
joins=c(length(x1))



powerpost(40,0.99,0.01,iterations,r,f,y,4000,joins)