###############################################################################################
###############################################################################################
#This function is optimized for the case when all intervalar lower/upper censoring limits
#are finite
###############################################################################################
###############################################################################################

Kan.IC = function(a,b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  a1 = a-mu
  b1 = b-mu
  p = pmvnorm(lower = a,upper = b,mean = mu,sigma = Sigma)
  run = qfun(a = a1,b = b1,Sigma = Sigma)
  qa = run$qa
  qb = run$qb
  q = qa-qb
  muY = mu+ Sigma%*%q/p
  D = matrix(0,n,n)
  for(i in seqq){
    D[i,i] = a[i]*qa[i]
    D[i,i] = D[i,i]-b[i]*qb[i]
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    ma = mu[-i]+Sigma[-i,i]/Sigma[i,i]*a1[i]
    run1 = qfun(a[-i]-ma,b[-i]-ma,RR)
    qa1 = run1$qa
    qb1 = run1$qb
    wa = qa[i]*ma+dnorm(x = a[i],mean = mu[i],sd = s[i])*RR%*%(qa1-qb1)
    mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
    run2 = qfun(a[-i]-mb,b[-i]-mb,RR)
    qa2 = run2$qa
    qb2 = run2$qb
    wb = qb[i]*mb + dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%(qa2-qb2)
    D[i,-i] = wa-wb
  }
  varY = Sigma + Sigma%*%(D - q%*%t(muY))/p
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)
  return(list(mean = round(muY,4),EYY = round(EYY,4),varcov = round(varY,4)))
}

###############################################################################################
###############################################################################################
#This function is optimized for the case when it DOES exist infinite values in the lower/upper
#truncation limits
###############################################################################################
###############################################################################################

Kan.LRIC = function(a,b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  a1 = a-mu
  b1 = b-mu
  p = pmvnorm(lower = a,upper = b,mean = mu,sigma = Sigma)
  run = qfun(a = a1,b = b1,Sigma = Sigma)
  qa = run$qa
  qb = run$qb
  q = qa-qb
  muY = mu+ Sigma%*%q/p
  D = matrix(0,n,n)
  for(i in seqq){
    if(a[i] != -Inf){
      D[i,i] = a[i]*qa[i]
    }
    if(b[i] != Inf){
      D[i,i] = D[i,i]-b[i]*qb[i]
    }
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    if(a[i] == -Inf){
      wa = matrix(0,n-1,1)
    }else
    {
      ma = mu[-i]+Sigma[-i,i]/Sigma[i,i]*a1[i]
      run1 = qfun(a[-i]-ma,b[-i]-ma,RR)
      qa1 = run1$qa
      qb1 = run1$qb
      wa = qa[i]*ma+dnorm(x = a[i],mean = mu[i],sd = s[i])*RR%*%(qa1-qb1)
    }
    if(b[i] == Inf){
      wb = matrix(0,n-1,1)
    }else
    {
      mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
      run2 = qfun(a[-i]-mb,b[-i]-mb,RR)
      qa2 = run2$qa
      qb2 = run2$qb
      wb = qb[i]*mb + dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%(qa2-qb2)
    }
    D[i,-i] = wa-wb
  }
  varY = Sigma + Sigma%*%(D - q%*%t(muY))/p
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)
  return(list(mean = round(muY,4),EYY = round(EYY,4),varcov = round(varY,4)))
}

###############################################################################################
###############################################################################################
#This is the original Vaida's function for upper truncation (right censoring)
###############################################################################################
###############################################################################################

Kan.RC = function(b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  b1 = b-mu
  p = pmvnorm(upper = as.numeric(b),mean = as.numeric(mu),sigma = Sigma)
  qb = qfun_b(b1 = b1,Sigma = Sigma)
  muY = mu - Sigma%*%qb/p
  D = matrix(0,n,n)
  for(i in seqq){
    D[i,i] = D[i,i]-b[i]*qb[i]
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
    qb2 = qfun_b(b[-i]-mb,RR)
    wb = qb[i]*mb - dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%qb2
    D[i,-i] = -wb
  }
  varY = Sigma + Sigma%*%(D + qb%*%t(muY))/p
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)
  return(list(mean = round(muY,4),EYY = round(EYY,4),varcov = round(varY,4)))
}

# ########################
# #TESTING
# ########################
# 
# p = 4
# mu  = c(matrix(c(1:p/10)))
# s  = 2*matrix(rnorm(p^2),p,p)
# Sigma = S = round(s%*%t(s)/10,1)
# lambda = seq(-1,2,length.out = p)
# tau = 1
# 
# a = rep(-Inf,p)
# b = mu+2
# 
# Kan.R(b,mu,Sigma)
# Vaida(b,mu,Sigma)
# 
# compare <- microbenchmark(Kan.R(b,mu,Sigma),
#                           Vaida(b,mu,Sigma),
#                           times = 100)
# autoplot(compare)
# 
# ########################
# 
# a = mu-1
# b = mu+2
# a[2] = -Inf
# b[3] = Inf
# 
# Kan.LRIC(a,b,mu,Sigma)
# meanvarN(a,b,mu,Sigma)
