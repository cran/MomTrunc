onlymeanNuni = function(a=-Inf,b=Inf,mu,Sigma){
  p = length(mu)
  s = sqrt(Sigma)
  F0    = pnorm(b,mu,s)-pnorm(a,mu,s)
  da     = dnorm(a,mu,s)
  db     = dnorm(b,mu,s)
  F1     = mu*F0 + Sigma*(da - db)
  muY = c(F1/F0)
  return(list(mean = round(muY,4)))
}
