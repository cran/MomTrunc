#library(mnormt)
#library(mvtnorm)

#######################################################################################
#######################################################################################

sqrtm <- function(A)
{
  if(length(A)==1)
    Asqrt=sqrt(A)
  else{
    sva <- svd(A)
    if (min(sva$d)>=0){
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v)  # svd e decomposi??o espectral
    }else{
      stop("Matrix square root is not defined")
    }
  }
  return(as.matrix(Asqrt))
}

#######################################################################################
#######################################################################################

dmvESN1 <- function(y, mu=0, Sigma=1, lambda,tau){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- length(c(y))
  p <- 1
  tautil<-tau/sqrt(1+sum(lambda^2))
  dens <- dnorm(x = c(y),mean = c(mu),sd = sqrt(Sigma))*
    pnorm(apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
                  (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum)+tau)/
    pnorm(tautil)
  return(dens)
}

#######################################################################################
#######################################################################################

AcumESN<-function(y=c(1,1),mu=c(0,0),Sigma=diag(2),lambda=c(2,-1),tau=1){
  if(all(y == -Inf)){return(0)}
  if(all(y ==  Inf)){return(1)}
  y<-as.matrix(y)
  mu<-as.matrix(mu)
  lambda<-as.matrix(lambda)
  delta<-lambda/sqrt(1+sum(lambda^2))
  tautil<-tau/sqrt(1+sum(lambda^2))
  p<-length(y)
  yaum<- c(y-mu,tau)
  Omega1<- cbind(Sigma,-sqrtm(Sigma)%*%lambda)
  Omega2<- cbind(-t(sqrtm(Sigma)%*%lambda),1+sum(lambda^2))
  Omega<- rbind(Omega1,Omega2)
  rownames(Omega) <- colnames(Omega)
  resp<- pmvnorm(upper = yaum,sigma = Omega)/pnorm(tautil)
  return(resp[1])
}

#######################################################################################
#######################################################################################

pmvESN = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,tau){
  tautil<-tau/sqrt(1+sum(lambda^2))
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
  return(pmvnorm(lower = aaum,upper = baum,mean = rep(0,p+1),sigma = Omega)[1]/pnorm(tautil))
}

#######################################################################################
#######################################################################################

invmills = function(x,mu=0,sd=1){
  z = (x-mu)/sd
  if(z < -37.5){return(-z/sd)}else{return(dnorm(x,mu,sd)/pnorm(x,mu,sd))}
}
