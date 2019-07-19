#######################################################################################
#######################################################################################

rmvESN<-function(n,mu=c(0,0),Sigma=diag(2),lambda=c(-1,1),tau=1){
  #Validating Lambda
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector fot the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
    if(all(lambda==0)){
      warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
      out = rmvnorm(n = n,mean = mu,sigma = Sigma)
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the ESN case (zero for the Skew-normal case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
    tautil = tau/sqrt(1+sum(lambda^2))
    if(tautil < -2.2){
      if(tautil< -37){
        #print("normal aproximation")
        Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
        return(rmvnorm(n = n,mean = c(mu - tautil*Delta),sigma = Sigma - Delta%*%t(Delta)))
      }else{
        stop("Tau parameter too small leading to high rejection. Try incresing tau.")
      }
    }
    return(rESN0(n = n,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau))
  }
}

rESN0<-function(n = 10000,mu=c(0,2),Sigma=diag(2),lambda=c(-1,3),tau=1){
  mu<-as.matrix(mu)
  lambda<-as.matrix(lambda)
  varphi<-lambda/sqrt(1+sum(lambda^2))
  tautil<-tau/sqrt(1+sum(lambda^2))
  p<-length(mu)
  SS = sqrtm(Sigma)
  Omega1<- cbind(Sigma,-SS%*%varphi)
  Omega2<- cbind(-t(SS%*%varphi),1)
  Omega<- rbind(Omega1,Omega2)
  gen <- rmvnorm(n = round(n/max(pnorm(tautil)-0.01,0.02)),mean = rep(0,p+1),sigma = Omega)
  bool = gen[,p+1] < tautil
  nn = sum(bool)
  gen1 = gen[bool,1:p]
  muM = matrix(rep(mu,nn),nrow = nn,ncol = p,byrow = TRUE)
  gen2 = muM + gen1
  nout = min(n,nrow(gen2))
  return(gen2[1:nout,])
}
