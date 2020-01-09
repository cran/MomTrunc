cdfFN = function(x,mu,Sigma)
{
  p = length(mu)
  if(p==1)
  {
    return(pnorm2(-x,x,mu,sqrt(Sigma)))
  }else
  {
    return(pmvn.genz(lower = -x-mu,upper = x-mu,sigma = Sigma)$Estimation)
  }
}

