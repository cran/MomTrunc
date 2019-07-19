#library(mnormt)
#library(mvtnorm)

select = function(a,b,mu,s){
  a1 = (a-mu)/s
  b1 = (b-mu)/s
  if(a1<0 & b1<0){
    val = b - exp(b/10)
  }
  if(a1>0 & b1>0){
    val = a + exp(-a/10)
  }
  return(val)
}

#######################################################################################
#######################################################################################

pnorm2 = function(lower = -Inf,upper = Inf,mean = 0,sd = 1){
  if(lower == -Inf){
    return(pnorm(upper,mean = mean,sd = sd))
  }else{
    if(upper == Inf){
      return(pnorm(lower,mean = mean,sd = sd,lower.tail = FALSE))
    }else{
      return(pnorm(upper,mean = mean,sd = sd) - pnorm(lower,mean = mean,sd = sd))
    }
  }
}

#######################################################################################
#######################################################################################

is.pd = function(MAT){
  all(eigen(MAT)$values > 0)
}

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
  s = sqrt(Sigma)
  tautil = tau/sqrt(1+sum(lambda^2))
  if(tautil< -35){
    #print("normal aproximation")
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    return(dnorm(y,mu-mub,sqrt(Gamma)))
  }
  dens <- dnorm(x = c(y),mean = c(mu),sd = sqrt(Sigma))*
    pnorm(apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
                  (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum)+tau)/pnorm(tautil)
  return(dens)
}

#######################################################################################
#######################################################################################

AcumESN<-function(y=c(1,1),mu=c(0,0),Sigma=diag(2),lambda=c(2,-1),tau=1){
  if(all(y == -Inf)){return(0)}
  if(all(y ==  Inf)){return(1)}
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -35){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(pmvnorm(upper = y,mean = c(mu - tautil*Delta),sigma = Sigma - Delta%*%t(Delta))[1])
  }
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

invmills = function(x,mu=0,sd=1){
  z = (x-mu)/sd
  if(z < -37.5){return(-z/sd)}else{return(dnorm(x,mu,sd)/pnorm(x,mu,sd))}
}
