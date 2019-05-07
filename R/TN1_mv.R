meanvarNuni = function(a=-Inf,b=Inf,mu,Sigma){
  p = length(mu)
  s = sqrt(Sigma)
  F0    = pnorm(b,mu,s)-pnorm(a,mu,s)
  da     = dnorm(a,mu,s)
  db     = dnorm(b,mu,s)
  F1     = mu*F0 + Sigma*(da - db)
  muY = c(F1/F0)
  F2 = mu*F1 + Sigma*F0 + Sigma*(ifelse(a == -Inf,0,a*da) - ifelse(b == Inf,0,b*db))
  varY = c(F2/F0) - muY^2
  return(list(mean = round(muY,4),EYY = round(c(F2/F0),4), varcov = round(varY,4)))
}
