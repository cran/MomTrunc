#######################################################################################
#######################################################################################

pmvESN = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,tau, ...){
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -37){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    Gamma = Sigma - Delta%*%t(Delta)
    rownames(Gamma) <- colnames(Gamma)
    return(pmvnorm(lower = lower,upper = upper,mean = c(mu - tautil*Delta),sigma = Gamma, ...)[1])
  }
  aaum = c(lower-mu,-Inf)
  baum = c(upper-mu,tautil)
  mu<-as.matrix(mu)
  lambda<-as.matrix(lambda)
  varphi<-lambda/sqrt(1+sum(lambda^2))
  p<-length(mu)
  SS = sqrtm(Sigma)
  Omega1<- cbind(Sigma,-SS%*%varphi)
  Omega2<- cbind(-t(SS%*%varphi),1)
  Omega<- rbind(Omega1,Omega2)
  rownames(Omega) <- colnames(Omega)
  return(pmvnorm(lower = aaum,upper = baum,mean = rep(0,p+1),sigma = Omega, ...)[1]/pnorm(tautil))
}
