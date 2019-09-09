onlymeanN = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma){
  p = length(mu)
  if(p==1){
    out = onlymeanNuni(a = lower,b = upper,mu = mu,Sigma = Sigma)  #OK
    return(out)
  }
  if(all(is.infinite(c(lower,upper)))){
    return(list(mean = mu))
  }
  bool = is.infinite(lower) & is.infinite(upper)
  #if exists (-Inf,Inf) limits
  if(sum(bool)>0){
    out = withinfs_mean(lower,upper,mu,Sigma,bool)
  }else{
    out = Vaida.LRIC.onlymean(a = lower,b = upper,mu = mu,Sigma = Sigma)
  }
  return(out)
}
