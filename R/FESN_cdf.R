cdfFESN = function(x,mu,Sigma,lambda,tau)
{
  p = length(mu)
  if(p==1)
  {
    return(AcumESN(x,mu,Sigma,lambda,tau) - AcumESN(-x,mu,Sigma,lambda,tau))
  }else
  {
    return(pmvESN(lower = -x-mu,upper = x-mu,Sigma = Sigma,lambda = lambda,tau = tau))
  }
}
