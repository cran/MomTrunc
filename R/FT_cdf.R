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
