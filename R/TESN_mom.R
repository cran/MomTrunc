KmomentESN = function(k,a,b,mu,Sigma,lambda,tau)
{
    p = length(mu)
    SS   = sqrtm(Sigma)
    tautil = tau/sqrt(1+sum(lambda^2))
    varpsi = lambda/sqrt(1+sum(lambda^2))
    Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
    rownames(Omega) <- colnames(Omega)
    return(KmomentN(c(k,0),a = c(a,-10^7),b = c(b,tautil),mu = c(mu,0),Sigma = Omega)[,-(p+1)])
}
