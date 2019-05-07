###############################################################################################
###############################################################################################
#This function compute the mean and the var-cov matrix for a TN distribution
###############################################################################################
###############################################################################################

meanvarN7 = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma){
  p = length(mu)
  if(p==1){
    out = meanvarNuni(a = lower,b = upper,mu = mu,Sigma = Sigma)  #OK
    return(out)
  }
  if(all(is.infinite(lower))){
    if(all(is.infinite(upper))){
      #No truncating at all
      return(list(mean = mu,EYY = Sigma + mu%*%t(mu),varcov = Sigma))
    }else{
      #Right censoring
      if(p<4){
        out = Kan.RC(b = upper,mu = mu,Sigma = Sigma)
      }else{
        out = Vaida.RC(b = upper,mu = mu,Sigma = Sigma)   #OK
      }
    }
  }else{
    if(all(is.infinite(upper))){
      #Left censoring
      if(p<4){
        out = Kan.RC(b = -lower,mu = -mu,Sigma = Sigma) #OK
      }else{
        out = Vaida.RC(b = -lower,mu = -mu,Sigma = Sigma) #OK
      }
      out$mean = -out$mean
    }else{
      #intervalar censoring
      if(all(is.finite(c(lower,upper)))){
        #no infinites #all intervalar truncated
        if(p<4){
          out = Kan.IC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }else{
          out = Vaida.IC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }
      }else
      {
        #All kind of censoring
        if(p<4){
          out = Kan.LRIC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }else{
          out = Vaida.LRIC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }
      }
    }
  }
  return(out)
}

# ########################
# #TESTING
# ########################
#
# p = 3
# mu  = c(matrix(c(1:p/10)))
# s  = 2*matrix(rnorm(p^2),p,p)
# Sigma = S = round(s%*%t(s)/10,1)
# lambda = seq(-1,2,length.out = p)
# tau = 1
#
# ########################################################################
# #INTERVALAR TRUNCATION ALL FINITES LIMITS #OK
# ########################################################################
# a = mu-2
# #a = rep(-Inf,p)
# #b = rep(Inf,p)
# b = mu+1
# ########################
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
# meanvarN7(a,b,mu,Sigma)
#
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 5)
# autoplot(compare)
# 10.8/5.8
#
# ########################################################################
# #LEFT TRUNCATED #OK
# ########################################################################
# a = mu-2
# #a = rep(-Inf,p)
# b = rep(Inf,p)
# #b = mu+1
# ########################
# meanvarESN7(a,b,mu,Sigma,lambda,tau)
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 5)
# autoplot(compare)
# 2/1.18
#
# ########################################################################
# #RIGHT TRUNCATED #OK
# ########################################################################
# #a = mu-2
# a = rep(-Inf,p)
# #b = rep(Inf,p)
# b = mu+1
# ########################
# meanvarESN7(a,b,mu,Sigma,lambda,tau)
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 10)
# autoplot(compare)
# 2.01/1.14
#
# ########################################################################
# #NO TRUNCATION
# ########################################################################
# #a = mu-2
# a = rep(-Inf,p)
# b = rep(Inf,p)
# #b = mu+1
# ########################
# meanvarESN7(a,b,mu,Sigma,lambda,tau)
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 50)
# autoplot(compare)
# 1946/204
#
# ########################################################################
# #INTERVALAR TRUNCATION RANDOM TRUNCATION
# ########################################################################
# a = mu-2
# b = mu+1
# gens = round(runif(n = p,min = 1,max = p))
# sel1 = 1:round(p/2)
# sel1c = (round(p/2)+1):p
# a[gens[sel1]] = -Inf
# b[gens[sel1c]] = Inf
# ########################
# meanvarESN7(a,b,mu,Sigma,lambda,tau)
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 20)
# autoplot(compare)
# 330/190
# ########################################################################
#
# p = 50
# mu  = c(matrix(c(1:p/10)))
# s  = 2*matrix(rnorm(p^2),p,p)
# Sigma = S = round(s%*%t(s)/10,1)
# lambda = seq(-1,2,length.out = p)
# tau = 1
#
# a = mu-3
# b = mu+2
# gens = round(runif(n = p,min = 1,max = p))
# sel1 = 1:round(p/2)
# sel1c = (round(p/2)+1):p
# a[gens[sel1]] = -Inf
# b[gens[sel1c]] = Inf
# ########################
# meanvarESN7(a,b,mu,Sigma,lambda,tau)
# meanvarESN3(a,b,mu,Sigma,lambda,tau)
#
# compare <- microbenchmark(meanvarESN7(a,b,mu,Sigma,lambda,tau),
#                           meanvarESN3(a,b,mu,Sigma,lambda,tau),
#                           times = 2)
# autoplot(compare)
# 698/350
