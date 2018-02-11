# require(matrixcalc)
# require(mvtnorm)

####################################################################
#AUXILIAR CODES T
####################################################################

#location-scale student-T distribution pdf
dent<-function(x,mu,sigma2,nu){
  z<-(x-mu)/sqrt(sigma2)
  return(1/sqrt(sigma2)*dt(x = z,df = nu))
}

#########################################################################

#location-scale student-T distribution pcf
pent<-function(x,mu,sigma2,nu){
  return(pt((x-mu)/sqrt(sigma2),nu))
}

#########################################################################

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

#########################################################################

F1k = function(k,a,b,mu,sigma2,nu)
{
  if(k > nu-2) stop("Moment k must be grater than nu-2")
  if(k<0) return(0)
  if(k==0){return(list(ks = 0, ress = pent(b,mu,sigma2,nu) - pent(a,mu,sigma2,nu)))}

  nnu = nu/(nu-2)
  fk = F1k(k-1,a,b,mu,sigma2,nu)
  res = mu*fk$ress[1] + (k-1)*nnu*sigma2*F1kgen(k-2,a,b,mu,nnu*sigma2,nu-2) +
    nnu*sigma2*(a^(k-1)*dent(a,mu,nnu*sigma2,nu-2) - b^(k-1)*dent(b,mu,nnu*sigma2,nu-2))
  return(list(ks = c(k,fk$ks), ress = c(res,fk$ress)))
}

#########################################################################

F1kgen = function(k,a,b,mu,sigma2,nu)
{
  if(k > nu-2) stop("Moment k must be grater than nu-2")
  if(k<0) return(0)
  if(k==0){return(pent(b,mu,sigma2,nu) - pent(a,mu,sigma2,nu))}
  nnu = nu/(nu-2)
  res = mu*F1kgen(k-1,a,b,mu,sigma2,nu) + (k-1)*nnu*sigma2*F1kgen(k-2,a,b,mu,nnu*sigma2,nu-2) +
    nnu*sigma2*(a^(k-1)*dent(a,mu,nnu*sigma2,nu-2) - b^(k-1)*dent(b,mu,nnu*sigma2,nu-2))
  return(res)
}

#########################################################################

ckgen = function(p,eim,k,a,b,mu,Sigma,nu,cuts= 2.5e4)
{
  out = rep(NA,p)

  SSigma  = (nu/(nu-2))*Sigma
  ssigma2 = diag(SSigma)
  zA      = ((a - mu)^2)/diag(SSigma)
  zB      = ((b - mu)^2)/diag(SSigma)
  deltaA  = (nu - 2 + zA)/(nu - 1)
  deltaB  = (nu - 2 + zB)/(nu - 1)
  yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/diag(SSigma)

  for(j in 1:p)
  {
    mua   = mu[-j] + yA[j]*SSigma[,j][-j]
    mub   = mu[-j] + yB[j]*SSigma[,j][-j]
    SSS   = SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]

    out[j] = k[j]*Fkgen(k-eim[j,],a,b,mu,SSigma,nu-2) +
      a[j]^k[j]*dent(a[j],mu[j],ssigma2[j],nu-2)*Fkgen(k[-j],a[-j],b[-j],mua,deltaA[j]*SSS,nu-1,cuts) -
      b[j]^k[j]*dent(b[j],mu[j],ssigma2[j],nu-2)*Fkgen(k[-j],a[-j],b[-j],mub,deltaB[j]*SSS,nu-1,cuts)
  }
  return(out)
}

#########################################################################

Fkgen = function(k,a,b,mu,Sigma,nu,cuts= 2.5e4)
{
  p = length(k)
  GB = GenzBretz(maxpts = cuts,abseps = 1e-9,releps = 0)
  if(p==1)
  {
    return(as.numeric(F1kgen(k,a,b,mu,Sigma,nu)))
  }else{
    if(any(k<0)) return(0)
    if(all(k == 0)){return(as.numeric(pmvt(lower = a - mu,upper = b - mu,df = nu,sigma = Sigma, algorithm = GB)))}
    eim  = diag(p)
    nnu  = nu/(nu-2)
    i = min(seq(1,p)[k>0])
    kk = k-eim[i,]
    ck = ckgen(p,eim,kk,a,b,mu,Sigma,nu)
    res  = mu[i]*Fkgen(kk,a,b,mu,Sigma,nu) + as.numeric(nnu*as.matrix(t(eim[i,]))%*%Sigma%*%ck)
    return(res)
  }
}

#########################################################################

Fk = function(k,a,b,mu,Sigma,nu,cuts= 2.5e4)
{
  p = length(k)
  GB = GenzBretz(maxpts = cuts,abseps = 1e-9,releps = 0)
  if(p==1)
  {
    res = F1k(k,a,b,mu,Sigma,nu)
    return(list(ks = res$ks, ress = res$ress))
  }else{
    if(any(k<0)) return(list(ress= 0))
    if(all(k == 0)){return(list(ks = rep(0,p),ress = as.numeric(pmvt(lower = a - mu,upper = b - mu,df = nu,sigma = Sigma,algorithm = GB))))}
    eim  = diag(p)
    nnu  = nu/(nu-2)
    i = min(seq(1,p)[k>0])
    kk = k-eim[i,]
    ck = ckgen(p,eim,kk,a,b,mu,Sigma,nu,cuts)
    fk = Fk(kk,a,b,mu,Sigma,nu)
    res1 = fk$ress[1]
    res  = mu[i]*fk$ress[1] + as.numeric(nnu*as.matrix(t(eim[i,]))%*%Sigma%*%ck)
    return(list(ks = rbind(k,fk$ks), ress = c(res,fk$ress)))
  }
}

####################################################################
#AUXILIAR CODES NORMAL
####################################################################

qfun = function(a,b,Sigma)
{
  n = length(a)
  s = sqrt(diag(Sigma))
  if(n==1){
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa,qb = qb))
  }

  qa = qb = rep(0,n)
  for(i in 1:n){
    if(a[i] != -Inf){
      qa[i] = dnorm(x = a[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*a[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    }
    if(b[i] != Inf){
      qb[i] = dnorm(x = b[i],mean = 0,sd = s[i])*pmvnorm(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*b[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))
    }
  }
  return(list(qa = qa,qb = qb))
}


recintab = function(kappa,a,b,mu,Sigma)
{
  n = length(kappa)
  seqq = seq_len(n)
  if(n==1)
  {
    M = matrix(data = 0,nrow = kappa+1,ncol = 1)
    s1 = sqrt(Sigma)
    aa = (a-mu)/s1
    bb = (b-mu)/s1
    M[1] = pnorm(bb)-pnorm(aa)
    if(kappa>0)
    {
      pdfa = s1*dnorm(aa)
      pdfb = s1*dnorm(bb)
      M[2] = mu*M[1]+pdfa-pdfb
      if(a == -Inf){a = 0}
      if(b == Inf){b = 0}
      for(i in seq_len(kappa-1)+1)
      {
        pdfa = pdfa*a
        pdfb = pdfb*b
        M[i+1] = mu*M[i] + (i-1)*Sigma*M[i-1] + pdfa - pdfb
      }
    }
  }else
  {
    #
    #  Create a matrix M, with its nu-th element correpsonds to F_{nu-2}^n(mu,Sigma).
    #
    M  = array(data = 0,dim = kappa+1)
    pk = prod(kappa+1)
    #
    #  We create two long vectors G and H to store the two different sets
    #  of integrals with dimension n-1.
    #
    nn = round(pk/(kappa+1),0)
    begind = cumsum(c(0,nn))
    pk1 = begind[n+1]   # number of (n-1)-dimensional integrals
    #  Each row of cp corresponds to a vector that allows us to map the subscripts
    #  to the index in the long vectors G and H
    cp = matrix(0,n,n)
    for(i in seqq)
    {
      kk = kappa
      kk[i] = 0
      cp[i,] = c(1,cumprod(kk[1:n-1]+1))
    }
    G = rep(0,pk1)
    H = rep(0,pk1)
    s = sqrt(diag(Sigma))
    pdfa = dnorm(x = a,mean = mu,sd = s)
    pdfb = dnorm(x = b,mean = mu,sd = s)
    for(i in seqq)
    {
      ind2 = seqq[-i]
      kappai = kappa[ind2]
      ai = a[ind2]
      bi = b[ind2]
      mui = mu[ind2]
      Si = Sigma[ind2,i]
      SSi = Sigma[ind2,ind2] - Si%*%t(Si)/Sigma[i,i]
      ind = (begind[i]+1):begind[i+1]
      if(a[i] != -Inf){
        mai = mui + Si/Sigma[i,i]*(a[i]-mu[i])
        G[ind] = pdfa[i]*recintab(kappai,ai,bi,mai,SSi)
      }
      if(b[i] != Inf){
        mbi = mui+Si/Sigma[i,i]*(b[i]-mu[i])
        H[ind] = pdfb[i]*recintab(kappai,ai,bi,mbi,SSi)
      }
    }
    M[1] = pmvnorm(lower=a, upper=b, mean=mu,sigma = Sigma)
    a[a == -Inf] = 0
    b[b == Inf] = 0
    cp1 = cp[n,]
    dims  = lapply(X = kappa + 1,FUN = seq_len)
    grid = do.call(expand.grid, dims)

    for(i in seq_len(pk-1)+1)
    {
      kk = as.numeric(grid[i,])
      ii = sum((kk-1)*cp1)+1
      i1 = min(seqq[kk>1])
      kk1 = kk
      kk1[i1] = kk1[i1]-1
      ind3 = ii-cp1[i1]
      M[ii] = mu[i1]*M[ind3]
      for(j in seqq)
      {
        kk2 = kk1[j]-1
        if(kk2 > 0)
        {
          M[ii] = M[ii]+Sigma[i1,j]*kk2*M[ind3-cp1[j]]
        }
        ind4 = as.numeric(begind[j] + sum(cp[j,]*(kk1-1)) - cp[j,j]*kk2+1)
        M[ii] = M[ii]+Sigma[i1,j]*(a[j]^kk2*G[ind4]-b[j]^kk2*H[ind4])
      }
    }
  }
  return(M)
}



####################################################################
#AUXILIAR CODES FOLDED NORMAL
####################################################################

allFN = function(mu,Sigma){
  n = length(mu)
  a = rep(0,n)
  b = rep(Inf,n)
  s = sqrt(diag(Sigma))
  if(n==1){
    a1 = (a-mu)/s
    b1 = (b-mu)/s
    p = pnorm(b1)-pnorm(a1)
    muY = mu+(dnorm(a1)-dnorm(b1))*s/p
    if(a == -Inf){
      a = 0
    }
    if(b == Inf){
      b = 0
    }
    varY = Sigma+(mu-muY)*muY+(a*dnorm(a1)-b*dnorm(b1))*s/p
    return(list(F1 = p*muY,F2 = p^2*varY + muY^2))
  }
  seqq = seq_len(n)
  a1 = a-mu
  b1 = b-mu
  p = pmvnorm(lower = a,upper = b,mean = mu,sigma = Sigma)
  run = qfun(a = a1,b = b1,Sigma = Sigma)
  qa = run$qa
  qb = run$qb
  q = qa-qb
  muY = mu+ Sigma%*%q/p
  D = matrix(0,n,n)
  for(i in seqq){
    if(a[i] != -Inf){
      D[i,i] = a[i]*qa[i]
    }
    if(b[i] != Inf){
      D[i,i] = D[i,i]-b[i]*qb[i]
    }
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    if(a[i] == -Inf){
      wa = matrix(0,n-1,1)
    }else
    {
      ma = mu[-i]+Sigma[-i,i]/Sigma[i,i]*a1[i]
      run1 = qfun(a[-i]-ma,b[-i]-ma,RR)
      qa1 = run1$qa
      qb1 = run1$qb
      wa = qa[i]*ma+dnorm(x = a[i],mean = mu[i],sd = s[i])*RR%*%(qa1-qb1)
    }
    if(b[i] == Inf){
      wb = matrix(0,n-1,1)
    }else
    {
      mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
      run2 = qfun(a[-i]-mb,b[-i]-mb,RR)
      qa2 = run2$qa
      qb2 = run2$qb
      wb = qb[i]*mb + dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%(qa2-qb2)
    }
    D[i,-i] = wa-wb
  }
  varY = Sigma + Sigma%*%(D - q%*%t(muY))/p
  return(list(F1 = p*muY,F2 = p*(varY + muY%*%t(muY))))
}

####################################################################
#AUXILIAR CODES FOLDED STUDENT
####################################################################

allFT = function(mu,Sigma,nu)
{
  p = length(mu)
  nnu = nu/(nu-2)
  if(nu>=4){
    if(p==1){
      F0 = 1 - pent(0,mu,Sigma,nu)
      nnusigma2 = nnu*Sigma
      ta = dent(0,mu,nnusigma2,nu-2)
      F1 = mu*F0 + nnusigma2*(ta)
      F2 = mu*F1 + nnusigma2*(1 - pent(0,mu,nnusigma2,nu-2))
      return(list(F1 = F1,F2 = F2))
    }
    a = rep(0,p)
    b = rep(Inf,p)
    GB = GenzBretz(maxpts = 5e4, abseps = 1e-5, releps = 0)
    F0 = pmvt(lower = -mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB)[1]
    F0nnu = pmvt(lower = -mu,upper = b-mu,df = nu - 2,sigma = nnu*Sigma, algorithm = GB)[1]

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

      aux1     = onlymeanT(a = a[-j],b = b[-j],mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
      ca[j]    = dent(0,mu[j],ssigma2[j],nu-2)*aux1$F00
      Wa[-j,j] = aux1$muY
    }
    F11  = F0*mu + SSigma%*%ca
    F22  = F11%*%t(mu) +  (F0nnu*diag(p) + Wa%*%diag(ca))%*%SSigma
    F22 = (F22 + t(F22))/2
    return(list(F1 = F11,F2 = F22))
  }else{
    if(p==1){
      F0 = 1 - pent(0,mu,Sigma,nu)
      ta = dent(0,mu,nnu*Sigma,nu-2)
      F1 = mu*F0 + nnu*Sigma*(ta)
      return(list(F1 = F1))
    }
    a = rep(0,p)
    b = rep(Inf,p)
    GB = GenzBretz(maxpts = 5e4, abseps = 1e-5, releps = 0)
    F0 = pmvt(lower = -mu,upper = b-mu,df = nu,sigma = Sigma, algorithm = GB)[1]

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
     ca[j] = dent(0,mu[j],ssigma2[j],nu-2)*pmvt(lower = a[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),upper = b[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),df = nu-1,sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]), algorithm = GB)[1]
    }
    F11  = F0*mu + SSigma%*%ca
    return(list(F1 = F11))
  }
}
