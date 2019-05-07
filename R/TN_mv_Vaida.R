###############################################################################################
###############################################################################################
#This function is optimized for the case when all intervalar lower/upper censoring limits
#are finite
###############################################################################################
###############################################################################################

Vaida.IC = function(a=-c(3,2),b=c(1,2),mu=c(0,0),Sigma=matrix(c(1,-0.5,-0.5,1),2,2)) {
  p=length(mu)
  if(p==1){
    s = sqrt(Sigma)
    a1 = (a-mu)/s
    b1 = (b-mu)/s
    L = pnorm(b1)-pnorm(a1)
    muY = mu+(dnorm(a1)-dnorm(b1))*s/L
    varY = Sigma+(mu-muY)*muY+(a*dnorm(a1)-b*dnorm(b1))*s/L
    return(list(mean = round(muY,4),EYY = round(varY+muY^2,4),varcov = round(varY,4)))
  }else{
    a1 <- (a-mu)/sqrt(diag(Sigma))
    b1 <- (b-mu)/sqrt(diag(Sigma))
    R <-  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    alpha <- pmvnorm(lower = as.vector(a1),upper=as.vector(b1),corr=R)[1]
    qq = qfun.noinf(a1,b1,Sigma = R)
    da = qq$qa
    db = qq$qb
    H =  RH = matrix(0,p,p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(x = as.vector(a1),sigma=R) - dmvnorm(x = as.vector(c(a1[1],b1[2])),sigma=R) - dmvnorm(x = as.vector(c(b1[1],a1[2])),sigma=R) + dmvnorm(x = as.vector(b1),sigma=R)
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }else{
    for(s in 1:(p-1)){
        for(t in (s+1):p){
          #print(s);print(t)
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          Sj = R[-c(s,t), c(s,t), drop=F]%*%invR
          V  =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = dmvnorm(x = a1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
            pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%a1[c(s,t),drop=F]),sigma=V) - 
            dmvnorm(x = c(a1[s],b1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
            pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(a1[s],b1[t])),sigma=V) - 
            dmvnorm(x = c(b1[s],a1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
            pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(b1[s],a1[t])),sigma=V) +
            dmvnorm(x = b1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
            pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%b1[c(s,t),drop=F]),sigma=V)
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    a1[is.infinite(a1)] = 0
    b1[is.infinite(b1)] = 0
    h = (a1*da - b1*db) - apply(RH, 1, sum)
    diag(H) = h
    EX <- -(1/alpha)*R%*%(da - db)   # a vector with a length of p
    EXX <- R + 1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu
    varyic <- diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy <- varyic+Eycens%*%t(Eycens)
  }
  return(list(mean=round(Eycens,4),EYY=round(E2yy,4),varcov=round(varyic,4)))
}

###############################################################################################
###############################################################################################
#This function is optimized for the case when it DOES exist infinite values in the lower/upper
#truncation limits
###############################################################################################
###############################################################################################

Vaida.LRIC<-function(a=-c(-Inf,2),b=c(1,Inf),mu=c(0,0),Sigma=matrix(c(1,-0.5,-0.5,1),2,2)){
  p=length(mu)
  if(p==1){
    s = sqrt(Sigma)
    a1 = (a-mu)/s
    b1 = (b-mu)/s
    L = pnorm(b1)-pnorm(a1)
    muY = mu+(dnorm(a1)-dnorm(b1))*s/L
    if(a == -Inf){
      a = 0
    }
    if(b == Inf){
      b = 0
    }
    varY = Sigma+(mu-muY)*muY+(a*dnorm(a1)-b*dnorm(b1))*s/L
    return(list(mean = round(muY,4),EYY = round(varY+muY^2,4),varcov = round(varY,4)))
  }else{
    a1 <- (a-mu)/sqrt(diag(Sigma))
    b1 <- (b-mu)/sqrt(diag(Sigma))
    R <-  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    alpha <- pmvnorm(lower = as.vector(a1),upper=as.vector(b1),corr=R)[1]
    #mean
    qq = qfun(a1,b1,Sigma = R)
    da = qq$qa
    db = qq$qb
    #var
    H =  RH = matrix(0,p,p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(x = as.vector(a1),sigma=R) - dmvnorm(x = as.vector(c(a1[1],b1[2])),sigma=R) - dmvnorm(x = as.vector(c(b1[1],a1[2])),sigma=R) + dmvnorm(x = as.vector(b1),sigma=R)
      #sigma==R since b1 is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }else{
      for(s in 1:(p-1)){
        for(t in (s+1):p){
          #print(s);print(t)
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          Sj = R[-c(s,t), c(s,t), drop=F]%*%invR
          V  =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = ifelse(any(is.infinite(a1[c(s,t)])),0,dmvnorm(x = a1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
                                     pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%a1[c(s,t),drop=F]),sigma=V)) -
            ifelse(any(is.infinite(c(a1[s],b1[t]))),0,dmvnorm(x = c(a1[s],b1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
                     pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(a1[s],b1[t])),sigma=V)) -
            ifelse(any(is.infinite(c(b1[s],a1[t]))),0,dmvnorm(x = c(b1[s],a1[t]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
                     pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%c(b1[s],a1[t])),sigma=V)) +
            ifelse(any(is.infinite(b1[c(s,t)])),0,dmvnorm(x = b1[c(s,t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))*
                     pmvnorm(lower =as.vector(a1[-c(s,t)]),upper=as.vector(b1[-c(s,t)]),mean = as.vector(Sj%*%b1[c(s,t),drop=F]),sigma=V))
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    a1[is.infinite(a1)] = 0
    b1[is.infinite(b1)] = 0
    h = (a1*da - b1*db) - apply(RH, 1, sum)
    diag(H) = h
    EX <- -(1/alpha)*R%*%(da - db)   # a vector with a length of p
    EXX <- R + 1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu
    varyic <- diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy <- varyic+Eycens%*%t(Eycens)
  }
  return(list(mean=round(Eycens,4),EYY=round(E2yy,4),varcov=round(varyic,4)))
}


###############################################################################################
###############################################################################################
#This is the original Vaida's function for upper truncation (right censoring)
###############################################################################################
###############################################################################################

Vaida.RC = function(b=c(1,2),mu=c(0,0), Sigma = diag(2)){
  p=length(mu)
  if (p==1) {
    qq = (1/sqrt(Sigma))*(-b+mu)
    R = 1
    alpha = pnorm(-qq)
    dd = dnorm(-qq)
    H = qq*dd
    EX = (1/alpha)*dd   # a vector with a length of p
    EXX = 1+1/alpha*H
    varX = EXX-EX^2
    Eycens = -sqrt(Sigma)*EX+mu
    varyic= varX*Sigma
    E2yy=varyic+Eycens^2
  }
  else {
    qq = diag(1/sqrt(diag(Sigma)))%*%(-b+mu)
    R =  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    alpha = pmvnorm(upper=as.vector(-qq), corr=R)
    #print(qq)
    dd = rep(0, p)   #derivative vector
    for (j in 1:p){
      V = R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu = -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] = dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }
    H = matrix(rep(0, p*p), nrow=p)
    RH = matrix(rep(0, p*p), nrow=p)
    if(p==2){
      H[1,2] = H[2,1] = dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] = RH[2,1] = R[1,2]*H[1,2]
    }else{
      for( s in 1:(p-1)){
        for (t in (s+1):p){
          invR = solve(R[c(s,t), c(s,t), drop=F])
          nu = -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V =  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] = H[t,s] = pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] = RH[t,s] = R[s,t]*H[s,t]
        }
      }
    }
    h = qq*dd-apply(RH, 1, sum)
    diag(H) = h
    EX = (1/alpha)*R%*%dd
    EXX = R+1/alpha*R%*%H%*%R
    varX = EXX-EX%*%t(EX)
    Eycens = -diag(sqrt(diag(Sigma)))%*%EX+mu
    varyic = diag(sqrt(diag(Sigma)))%*%varX%*%diag(sqrt(diag(Sigma)))
    E2yy = varyic+Eycens%*%t(Eycens)
  }
  return(list(mean=round(Eycens,4),EYY=round(E2yy,4),varcov=round(varyic,4)))
}