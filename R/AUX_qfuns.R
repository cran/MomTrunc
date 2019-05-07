###################################################################################
###################################################################################
#AUXILIAR FUNCTIONS
###################################################################################
###################################################################################

qfun = function(a,b,Sigma)
{
  n = length(a)
  s = sqrt(diag(Sigma))
  if(n==1){
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa,qb = qb))
  }
  
  qa = qb = rep(0,n)
  for(i in 1:n){
    if(a[i] != -Inf){
      qa[i] = dnorm(x = a[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*a[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    }
    if(b[i] != Inf){
      qb[i] = dnorm(x = b[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*b[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    }
  }
  return(list(qa = qa,qb = qb))
}

###################################################################################
###################################################################################

qfun_b = function(b1,Sigma)
{
  p = length(b1)
  s = sqrt(diag(as.matrix(Sigma)))
  if(p==1){
    qb = dnorm(b1/s)/s
    return(qb)
  }
  qb = rep(0,p)
  for(i in 1:p){
    if(b1[i] != Inf){
      qb[i] = dnorm(x = b1[i],mean = 0,sd = s[i])*pmvnorm(upper = b1[-i],mean = Sigma[-i,i]/Sigma[i,i]*b1[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    }
  }
  return(qb)
}

###################################################################################
###################################################################################

qfun.noinf = function(a,b,Sigma)
{
  n = length(a)
  s = sqrt(diag(Sigma))
  if(n==1){
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa,qb = qb))
  }
  
  qa = qb = rep(0,n)
  for(i in 1:n){
    qa[i] = dnorm(x = a[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*a[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    qb[i] = dnorm(x = b[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*b[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
  }
  return(list(qa = qa,qb = qb))
}
