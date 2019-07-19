#########################################################################################################
#########################################################################################################

corrector = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bw=36){
  p = length(lower)
  ss = sqrt(diag(as.matrix(Sigma)))
  liminf = mu - bw*ss
  limsup = mu + bw*ss
  bool1 = upper < liminf
  bool2 = lower > limsup
  while(sum(bool1) +  sum(bool2) == 0){
    bw = bw - 2
    liminf = mu - bw*ss
    limsup = mu + bw*ss
    bool1 = upper < liminf
    bool2 = lower > limsup
  }
  bool   = bool1 | bool2
  seqq   = seq(1,p)
  bool   = seqq[bool]

  bool1  = seqq[bool1]
  val1   = upper[bool1]
  delta1 = exp(upper/10)
  val1   = val1 - delta1[bool1]

  bool2 = seqq[bool2]
  val2 = lower[bool2]
  delta2 = exp(-lower/10)
  val2 = val2 + delta2[bool2]

  val        = seqq
  val[bool1] = val1
  val[bool2] = val2
  val        = val[bool]

  if(length(bool) == p){
    out1 = matrix(val,p,1)
    out3 = matrix(0.0001,p,p)
    out2 = out3 + out1%*%t(out1)
    return(list(mean = out1,EYY = out2,varcov = out3))
  }
  tilSj = Sigma[-bool,-bool] - Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%Sigma[bool,-bool]
  op1 = meanvarN7(lower = lower[-bool],
                  upper = upper[-bool],
                  mu = c(mu[-bool] + Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%(val - mu[bool])),
                  Sigma = tilSj)
  op2 = op1
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = val
  op2$mean[-bool] = op1$mean
  op2$varcov = matrix(0.0001,p,p)
  op2$varcov[-bool,-bool] = op1$varcov
  op2$EYY = op2$varcov + op2$mean%*%t(op2$mean)
  return(op2)
}

#corrector(upper = c(0,-50,-50),mu=mu,Sigma=Sigma)

#########################################################################################################

withinfs = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bool){
  #bool = is.infinite(lower) & is.infinite(upper)
  p = length(bool)
  seqq   = seq(1,p)
  bool   = seqq[bool]
  q      = length(bool) #infinite dims
  if(q>1 & q<4){
    out.finite = Kan.LRIC(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool])
  }else{
    out.finite = Vaida.LRIC(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool])
  }
  op2 = list()
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  op2$EYY    = matrix(NA,p,p)
  op2$varcov = matrix(NA,p,p)
  op2$varcov[bool,bool]  = Sigma[bool,bool] - Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(diag(p-q) - out.finite$varcov%*%solve(Sigma[-bool,-bool]))%*%Sigma[-bool,bool]
  op2$varcov[-bool,-bool] = out.finite$varcov
  op2$varcov[bool,-bool] = Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%out.finite$varcov
  op2$varcov[-bool,bool] = t(op2$varcov[bool,-bool])
  op2$EYY = op2$varcov + op2$mean%*%t(op2$mean)
  return(op2)
}
