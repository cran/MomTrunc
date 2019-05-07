meanvarESNuni = function(a=-Inf,b=Inf,mu,Sigma,lambda,tau){
  p = length(mu)
    s = sqrt(Sigma)
    tautil = tau/sqrt(1+sum(lambda^2))
    phi    = lambda/sqrt(1+sum(lambda^2))
    eta    = invmills(tau,0,sqrt(1+lambda^2))
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    F0N    = pnorm(b,mu-mub,sqrt(Gamma))-pnorm(a,mu-mub,sqrt(Gamma))
    F0     = AcumESN(b,mu,Sigma,lambda,tau) - AcumESN(a,mu,Sigma,lambda,tau)
    da     = dmvESN1(a,mu,Sigma,lambda,tau)
    db     = dmvESN1(b,mu,Sigma,lambda,tau)
    F1     = mu*F0 + Sigma*(da - db) + lambda*s*eta*F0N
    muY = c(F1/F0)
    #delta = lambda*s*eta
    F1N = (mu-mub)*F0N + Gamma*(dnorm(a,mu-mub,sqrt(Gamma)) - dnorm(b,mu-mub,sqrt(Gamma)))
    F2 = mu*F1 + Sigma*F0 + Sigma*(ifelse(a == -Inf,0,a*da) - ifelse(b == Inf,0,b*db)) + lambda*s*eta*F1N
    varY = c(F2/F0) - muY^2
    return(list(mean = round(muY,4),EYY = round(c(F2/F0),4), varcov = round(varY,4)))
}
