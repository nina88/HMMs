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

new_label_switch=function(RT, s, r, n, lambda, P, f)
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




