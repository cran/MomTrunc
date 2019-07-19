####################################################################
#AUXILIAR CODES T
####################################################################

onlymeanT = function(a,b,mu,Sigma,nu)
{
  p = length(mu)
  nnu = nu/(nu-2)
  if(p==1){
    F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
    nnusigma2 = nnu*Sigma
    ta = dent(a,mu,nnusigma2,nu-2)
    tb = dent(b,mu,nnusigma2,nu-2)
    F1 = mu*F0 + nnusigma2*(ta-tb)
    return(list(muY = as.numeric(F1/F0),F00 = F0))
  }
  if(p==2){
    GB = GenzBretz(maxpts = 2e4,abseps = 1e-5,releps = 0)
    F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB)[1]

    #Vectors ca and cb
    SSigma  = nnu*Sigma
    ssigma2 = diag(SSigma)
    ca = cb = rep(0,p)
    deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
    deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
    yA      = (a - mu)/diag(SSigma)
    yB      = (b - mu)/diag(SSigma)

    for(j in 1:p)
    {
      if(a[j]!=-Inf){
        ca[j] = dent(a[j],mu[j],ssigma2[j],nu-2)*(pent(b[-j],mu[-j] + yA[j]*SSigma[,j][-j],deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu-1) - pent(a[-j],mu[-j] + yA[j]*SSigma[,j][-j],deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu-1))
      }
      if(b[j]!= Inf){
        cb[j] = dent(b[j],mu[j],ssigma2[j],nu-2)*(pent(b[-j],mu[-j] + yB[j]*SSigma[,j][-j],deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu-1) - pent(a[-j],mu[-j] + yB[j]*SSigma[,j][-j],deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu-1))
      }
    }
  }
  if(p>2){
    #if(p<5){mp = 30000}else{mp = 50000}
    GB0 = GenzBretz(maxpts = p*1e4,abseps = 1e-5,releps = 0)
    F0  = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB0)[1]
    #GB1 = GenzBretz(maxpts = 5e4,abseps = 5e-5,releps = 0)

    #Vectors ca and cb
    SSigma  = nnu*Sigma
    ssigma2 = diag(SSigma)
    ca = cb = rep(0,p)
    deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
    deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
    yA      = (a - mu)/diag(SSigma)
    yB      = (b - mu)/diag(SSigma)

    for(j in 1:p)
    {
      if(a[j]!=-Inf){
        ca[j] = dent(a[j],mu[j],ssigma2[j],nu-2)*pmvt(lower = a[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),upper = b[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),df = nu - 1,sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]), algorithm = GB0)[1]
      }
      if(b[j]!= Inf){
        cb[j] = dent(b[j],mu[j],ssigma2[j],nu-2)*pmvt(lower = a[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),upper = b[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),df = nu - 1,sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]), algorithm = GB0)[1]
      }
    }
  }
  muY = mu + SSigma%*%(ca - cb)/F0
  return(list(muY = muY,F00 = F0))
}
