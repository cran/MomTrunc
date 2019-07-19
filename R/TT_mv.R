#######################################################################
#STUDENT T
#######################################################################

meanvarT = function(a,b,mu,Sigma,nu)
{
  if(nu>=4){
    p = length(mu)
    nnu = nu/(nu-2)
    if(p==1){
      F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
      nnusigma2 = nnu*Sigma
      ta = dent(a,mu,nnusigma2,nu-2)
      tb = dent(b,mu,nnusigma2,nu-2)
      F1 = mu*F0 + nnusigma2*(ta-tb)
      F2 = mu*F1 + nnusigma2*(pent(b,mu,nnusigma2,nu-2) - pent(a,mu,nnusigma2,nu-2) + ifelse(a==-Inf,0,a*ta) - ifelse(b==Inf,0,b*tb))
      return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
    }
    GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
    #print(GB$maxpts)
    F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB)[1]
    F0nnu = pmvt(lower = a-mu,upper = b-mu,df = nu - 2,sigma = nnu*Sigma, algorithm = GB)[1]

    #Vectors ca and cb
    SSigma  = nnu*Sigma
    ssigma2 = diag(SSigma)
    ca = cb = a0 = a1 = rep(0,p)
    deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
    deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
    yA      = (a - mu)/diag(SSigma)
    yB      = (b - mu)/diag(SSigma)

    Wa = Wb = matrix(0,p,p)

    for(j in 1:p)
    {
      #W matrix construction

      if(a[j]!=-Inf){
        aux1     = onlymeanT(a = a[-j],b = b[-j],mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
        ca[j]    = dent(a[j],mu[j],ssigma2[j],nu-2)*aux1$F00
        Wa[-j,j] = aux1$muY
        Wa[j,j]  = a[j]
      }
      if(b[j]!= Inf){
        aux2     = onlymeanT(a = a[-j],b = b[-j],mu = mu[-j] + yB[j]*SSigma[,j][-j],Sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
        cb[j]    = dent(b[j],mu[j],ssigma2[j],nu-2)*aux2$F00
        Wb[-j,j] = aux2$muY
        Wb[j,j]  = b[j]
      }
    }
    muY  = mu + SSigma%*%(ca - cb)/F0
    Exx  = muY%*%t(mu) +  (F0nnu*diag(p) + Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma/F0
    varY = Exx - muY%*%t(muY)
    varY = (varY + t(varY))/2
    return(list(mean = muY,EYY = Exx,varcov = varY))
  }else{
    p = length(mu)
    nnu = nu/(nu-2)
    if(p==1){
      F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
      ta = dent(a,mu,nnu*Sigma,nu-2)
      tb = dent(b,mu,nnu*Sigma,nu-2)
      F1 = mu*F0 + nnu*Sigma*(ta-tb)
      return(list(mean = round(F1/F0,4),EYY = NA,varcov = NA))
    }
    GB = GenzBretz(maxpts = p*1e4, abseps = 1e-6, releps = 0)
    #print(GB$maxpts)
    F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB)[1]

    #Vectors ca and cb
    SSigma  = nnu*Sigma
    ssigma2 = diag(SSigma)
    ca = cb = a0 = a1 = rep(0,p)
    deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
    deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
    yA      = (a - mu)/diag(SSigma)
    yB      = (b - mu)/diag(SSigma)

    for(j in 1:p)
    {
      if(a[j]!=-Inf){
        ca[j] = dent(a[j],mu[j],ssigma2[j],nu-2)*pmvt(lower = a[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),upper = b[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),df = nu - 1,sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]), algorithm = GB)[1]
      }
      if(b[j]!= Inf){
        cb[j] = dent(b[j],mu[j],ssigma2[j],nu-2)*pmvt(lower = a[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),upper = b[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),df = nu - 1,sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]), algorithm = GB)[1]
      }
    }
    muY  = mu + SSigma%*%(ca - cb)/F0
    return(list(mean = muY,EYY = matrix(NA,p,p),varcov = matrix(NA,p,p)))
  }
}
