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
      return(list(mean = round(F1/F0,4),varcov = round(F2/F0 - (F1/F0)^2,4)))
    }
    pb <- txtProgressBar(min = 1,max = p,style = 3)
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
      setTxtProgressBar(pb, j)
    }
    muY  = mu + SSigma%*%(ca - cb)/F0
    Exx  = muY%*%t(mu) +  (F0nnu*diag(p) + Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma/F0
    varY = Exx - muY%*%t(muY)
    varY = (varY + t(varY))/2
    close(pb)
    return(list(mean = round(muY,4),varcov = round(varY,4)))
  }else{
    p = length(mu)
    nnu = nu/(nu-2)
    if(p==1){
      F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
      ta = dent(a,mu,nnu*Sigma,nu-2)
      tb = dent(b,mu,nnu*Sigma,nu-2)
      F1 = mu*F0 + nnu*Sigma*(ta-tb)
      return(list(mean = round(F1/F0,4)))
    }
    pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
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
      setTxtProgressBar(pb, j)
    }
    muY  = mu + SSigma%*%(ca - cb)/F0
    close(pb)
    return(list(mean = round(muY,4)))
  }
}


KmomentT = function(k,a,b,mu,Sigma,nu,cuts = 2.5e4)
{
  p = length(k)
  if(all(k == 0))
  {
    res  = pmvt(lower = a - mu,upper = b - mu,df = nu,sigma = Sigma)[1]
    mat  = data.frame(cbind(rep(0,p),round(res,4),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)
  }else{
    res  = Fk(k,a,b,mu,Sigma,nu,cuts)
    mat  = cbind(res$ks,c(res$ress))
    F0   = mat[dim(mat)[1],dim(mat)[2]]
    mat  = data.frame(cbind(round(mat,4),round(mat[,dim(mat)[2]]/F0,4)),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)}
}

#######################################################################
#NORMAL
#######################################################################

meanvarN = function(a,b,mu,Sigma){
  n = length(mu)
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
    return(list(mean = round(muY,4),varcov = round(varY,4)))
  }
  pb <- txtProgressBar(min = 1,max = n,style = 3)
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
    setTxtProgressBar(pb, i)
  }
  varY = Sigma + Sigma%*%(D - q%*%t(muY))/p
  varY = (varY + t(varY))/2
  close(pb)
  return(list(mean = round(muY,4),varcov = round(varY,4)))
}

KmomentN = function(k,a,b,mu,Sigma)
{
  p = length(k)
  if(p==1){
    res  = rev(recintab(k,a,b,mu,Sigma))
    res2 = res/as.numeric(res)[1]
    mat  = data.frame(cbind(seq(k,0),res,res2),row.names = NULL)
    colnames(mat) = c("k","F(k)","E[k]")
    return(mat)
  }
  else{
    t = rep(NA,p)
    for(i in 1:p){
      t[i] = paste("k",i,"=c(",paste(seq(0,k[i]),collapse = ","),")",sep="",collapse = ",")
    }
    res  = recintab(k,a,b,mu,Sigma)
    res2 = res/as.numeric(res)[1]
    eval(parse(text = paste("dimnames(res) = dimnames(res2) = list(",paste(t,collapse = ","),")",sep = "")))

    df1 = as.data.frame.table(round(res,4))
    df2 = as.data.frame.table(round(res2,4))
    mat1 = data.matrix(df1[prod(k+1):1,], rownames.force = NA)
    mat2 = data.matrix(df2[prod(k+1):1,-(1:p)], rownames.force = NA)
    mat = cbind(mat1,mat2)
    mat[,1:p] = mat[,1:p] - 1
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)
  }
}
#######################################################################
#FOLDED NORMAL
#######################################################################

meanvarFN = function(mu,Sigma)
{
  p = length(mu)
  if(p==1)
  {
    muss = mu/sqrt(Sigma)
    muFT  = c(sqrt(Sigma)*(muss*(2*pnorm(muss)-1) + 2*dnorm(muss)))
    exxFT = c(Sigma + mu^2)
    varFT = exxFT - muFT^2
    return(list(mean = muFT,varcov = varFT))
  }
  pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
  signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
  means = matrix(NA,p,2^p)
  exxs  = vars = array(NA,dim = c(p,p,2^p))

  for(i in 1:(2^p))
  {
    Lambda   = diag(signs[i,])
    mus      = c(Lambda%*%mu)
    Ss       = Lambda%*%Sigma%*%Lambda
    res      = allFN(mus,Ss)
    means[,i] = res$F1
    exxs[,,i] = res$F2
    setTxtProgressBar(pb, i)
  }
  muFN  = apply(X = means,MARGIN = 1,FUN = sum)
  exxFN = matrix(apply(X = exxs,MARGIN = 1:2,FUN = sum),p,p)
  varFN = exxFN - muFN%*%t(muFN)
  close(pb)
  return(list(mean = round(muFN,4),varcov = round(varFN,4)))
}

KmomentFN = function(k,mu,Sigma)
{
  p = length(mu)
  if(all(k == 0)){
    mat  = data.frame(cbind(t(rep(0,p)),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }
  if(p==1)
  {
    mat  = data.frame(cbind(seq(k,0),as.matrix(KmomentN(k,a = 0,b = Inf,-mu,Sigma))[,2]+as.matrix(KmomentN(k,a = 0,b = Inf,mu,Sigma))[,2]),row.names = NULL)
    colnames(mat) = c("k","E[k]")
    return(mat)
  }else{
    pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
    signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
    means = matrix(NA,p,2^p)
    ressum  = array(NA,dim = c(prod(k+1),2^p))

    for(i in 1:(2^p))
    {
      Lambda   = diag(signs[i,])
      mus      = c(Lambda%*%mu)
      Ss       = Lambda%*%Sigma%*%Lambda
      mom       = as.matrix(KmomentN(k,a = rep(0,p),b = rep(Inf,p),mus,Ss))[,-(p+2)]
      ressum[,i] = as.matrix(mom[,-(1:p)])
      setTxtProgressBar(pb, i)
    }

    mat  = data.frame(cbind(mom[,1:p],matrix(apply(X = ressum,MARGIN = 1,FUN = sum),prod(k+1),1)),row.names = NULL)
    mat[prod(k+1),p+1] = 1
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    close(pb)
    return(mat)
  }
}

#######################################################################
#FOLDED T
#######################################################################

meanvarFT = function(mu,Sigma,nu)
{
  p = length(mu)
  if(nu>=4){
    if(p==1)
    {
      nnu = nu/(nu-2)
      nnusigma2 = nnu*Sigma
      muFT  = mu*(1 - 2*pent(0,mu,Sigma,nu)) + 2*nnusigma2*dent(0,mu,nnusigma2,nu-2)
      exxFT = mu^2 + nnusigma2
      varFT = exxFT - muFT^2
      #nnu2 = (nu-2)/(nu-4)
      #3th momment -> mu^2*muFT + 3*mu*nnusigma2*(1 - 2*pent(0,mu,nnusigma2,nu-2)) + 4*nnu2*nnusigma2^2*dent(0,mu,nnu2*nnusigma2,nu-4)
      #4th momment -> mu^4 + 6*mu^2*nnusigma2 + 3*nnu2*nnusigma2^2
      return(list(mean = muFT,varcov = varFT))
    }else
    {
      pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
      signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
      means = matrix(NA,p,2^p)
      exxs  = vars = array(NA,dim = c(p,p,2^p))

      for(i in 1:(2^p))
      {
        Lambda   = diag(signs[i,])
        mus      = c(Lambda%*%mu)
        Ss       = Lambda%*%Sigma%*%Lambda
        res      = allFT(mus,Ss,nu)
        means[,i] = res$F1
        exxs[,,i] = res$F2
        setTxtProgressBar(pb, i)
      }
    }
    muFT  = apply(X = means,MARGIN = 1,FUN = sum)
    exxFT = matrix(apply(X = exxs,MARGIN = 1:2,FUN = sum),p,p)
    varFT = exxFT - muFT%*%t(muFT)

    #third = mu^2*(muFT) +
    #  3*mu*nnusigma2*(1 - 2*pent(0,mu,nnusigma2,nu-2)) + 4*nnu2*nnusigma2^2*dent(0,mu,nnu2*nnusigma2,nu-4)
    close(pb)
    return(list(mean = round(muFT,4),varcov = round(varFT,4)))
  }else{
    if(p==1)
    {
      nnu = nu/(nu-2)
      nnusigma2 = nnu*Sigma
      muFT  = mu*(1 - 2*pent(0,mu,Sigma,nu)) + 2*nnusigma2*dent(0,mu,nnusigma2,nu-2)
      return(list(mean = round(muFT,4)))
    }else
    {
      pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
      signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
      means = matrix(NA,p,2^p)

      for(i in 1:(2^p))
      {
        Lambda   = diag(signs[i,])
        mus      = c(Lambda%*%mu)
        Ss       = Lambda%*%Sigma%*%Lambda
        res      = allFT(mus,Ss,nu)
        means[,i] = res$F1
        setTxtProgressBar(pb, i)
      }
    }
    muFT  = apply(X = means,MARGIN = 1,FUN = sum)
    close(pb)
    return(list(mean = round(muFT,4)))
  }
}

KmomentFT = function(k,mu,Sigma,nu)
{
  p = length(mu)
  if(all(k == 0)){
    mat  = data.frame(cbind(t(rep(0,p)),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }
  if(p==1)
  {
    mat  = data.frame(cbind(seq(k,0),as.matrix(KmomentT(k,a = 0,b = 10^10,-mu,Sigma,nu))[,2]+as.matrix(KmomentT(k,a = 0,b = 10^10,mu,Sigma,nu))[,2]),row.names = NULL)
    colnames(mat) = c("k","E[k]")
    return(mat)

  }else{
    pb <- txtProgressBar(min = 1,max = 2^p,style = 3)
    signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
    means = matrix(NA,p,2^p)
    ressum  = array(NA,dim = c(sum(k)+1,2^p))
    for(i in 1:(2^p))
    {
      Lambda   = diag(signs[i,])
      mus      = c(Lambda%*%mu)
      Ss       = Lambda%*%Sigma%*%Lambda
      mom       = as.matrix(KmomentT(k,a = rep(0,p),b = rep(10^10,p),mus,Ss,nu))[,-(p+2)]
      ressum[,i] = as.matrix(mom[,-(1:p)])
      setTxtProgressBar(pb, i)
    }
    mat  = data.frame(cbind(mom[,1:p],matrix(apply(X = ressum,MARGIN = 1,FUN = sum),sum(k)+1,1)),row.names = NULL)
    mat[sum(k)+1,p+1] = 1
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    close(pb)
    return(mat)
  }
}
############################################################################

cdfFT = function(x,mu,Sigma,nu)
{
  p = length(mu)
  if(p==1)
  {
    return(pent(x,mu,Sigma,nu) - pent(-x,mu,Sigma,nu))
  }else
  {
    return(pmvt(lower = -x-mu,upper = x-mu,df = nu,sigma = Sigma)[1])
  }
}
