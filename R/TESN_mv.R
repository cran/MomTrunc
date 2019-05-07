###############################################################################################
###############################################################################################
#This function compute the mean and the var-cov matrix for a TESN distribution using the
#Normal augmented model
###############################################################################################
###############################################################################################

meanvarESN7 = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma,lambda,tau){
  p = length(mu)
  if(p==1){
    out = meanvarESNuni(a = lower,b = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)  #OK
    return(out)
  }
  if(all(is.infinite(lower))){
    if(all(is.infinite(upper))){
      #No truncating at all
      return(ESN.NOTRUNC(mu,Sigma,lambda,tau))  #OK
    }else{
      #Right censoring
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)
      if(p<4){
        out = Kan.RC(b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
      }else{
        out = Vaida.RC(b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)   #OK
      }
    }
  }else{
    if(all(is.infinite(upper))){
      #Left censoring
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,t(varpsi)%*%SS),rbind(SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)
      if(p<4){
        out = Kan.RC(b = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega) #OK
      }else{
        out = Vaida.RC(b = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega) #OK
      }
      out$mean = -out$mean
    }else{
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)
      #intervalar censoring
      if(all(is.finite(c(lower,upper)))){
        #no infinites #all intervalar truncated
        if(p<4){
          out = Kan.IC(a = c(lower,-10^7),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }else{
          out = Vaida.IC(a = c(lower,-10^7),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }
      }else
      {
        #All kind of censoring
        if(p<4){
          out = Kan.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }else{
          out = Vaida.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }
      }
    }
  }
  return(list(mean = matrix(out$mean[-(p+1)],p,1),EYY = out$EYY[-(p+1),-(p+1)],varcov = out$varcov[-(p+1),-(p+1)]))
}
