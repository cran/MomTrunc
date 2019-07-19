#######################################################################################
#######################################################################################

dmvESN<-function(x,mu=c(0,0),Sigma=diag(2),lambda=c(-1,1),tau=1){
  #Validating Lambda
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector fot the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda)){
      stop("Lambda must be numeric and have same dimension than mu.")
    }
    if(all(lambda==0)){
      warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
      out = dmvnorm(x = x,mean = mu,sigma = Sigma)
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the ESN case (zero for the Skew-normal case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1){
      stop("Tau must be numeric real number.")
    }
    tautil = tau/sqrt(1+sum(lambda^2))
    if(tautil< -37){
      #print("normal aproximation")
      Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
      return(dmvnorm(x = x,mean = c(mu - tautil*Delta),sigma = Sigma - Delta%*%t(Delta)))
    }
    return(dmvESN0(x,mu,Sigma,lambda,tau))
  }
}

dmvESN0 <- function(y, mu, Sigma, lambda,tau){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  if(length(c(mu)) == 1){return(dmvESN1(c(y),mu,Sigma,lambda,tau))}
  if(is.matrix(y)){
    y = matrix(y,ncol = nrow(Sigma),byrow = TRUE)
  }
  n <- nrow(y)
  p <- ncol(y)
  tautil<-tau/sqrt(1+sum(lambda^2))
  dens <- dmvnorm(y, c(mu),Sigma)
  #*
   # pnorm(apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
    #              (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum)+tau)/pnorm(tautil)
  return(dens)
}