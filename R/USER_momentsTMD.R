#FULL MOMENTS TRUNCATED

momentsTMD = function(kappa,lower = NULL,upper = NULL,mu,Sigma,lambda = NULL, tau = NULL,dist,nu = NULL)
{
  mu = c(mu)
  lambda = c(lambda)
  #Validating dims data set
  if(ncol(as.matrix(kappa)) > 1 | !all(kappa >= 0) | length(c(kappa)) != length(c(mu))) stop("kappa must be numeric with same dimensions than mu.")
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("mu must be numeric and have just one column")

  #validate mean an Sigma dimensions

  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.pd(Sigma))stop("Sigma must be a square symmetrical real posite definite matrix.")
  }
  if(all(is.null(lower))){
    lower = rep(-Inf,length(mu))
  }else{
    if(length(c(lower)) != length(c(mu)) | !is.numeric(lower))stop("Lower bound must be numeric and have same dimension than mu.")
  }
  if(all(is.null(upper))){
    upper = rep(Inf,length(mu))
  }else{
    if(length(c(upper)) != length(c(mu)) | !is.numeric(upper))stop("Upper bound must be numeric and have same dimension than mu.")
  }
  if(all(lower < upper) == FALSE)stop("Lower bound must be lower than or equal to upper bound.")

  #validating distributions and nu parameter


  if(dist=="normal"){
    out = KmomentN(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma)
  }else
  {
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else
      {
        if(nu%%1!=0){
          stop("Degrees of freedom 'nu' must be an integer greater than 2.")
        }else
        {
          if(!all(kappa == 0) & nu < max(3,sum(kappa)+2)){stop("Sorry, we can only compute the kappa-th moment for degrees of freedom greater than or equal to 'sum(kappa)+2'.")
          }else
          {
            if(nu >= 200){
              warning("For degrees of freedom >= 200, Normal case is considered.",immediate. = TRUE)
              out = KmomentN(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma)
            }else
            {
              out = KmomentT(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma,nu = nu)
            }
          }
        }
      }
    }else
    {
      if(dist == "ESN" | dist == "SN"){
        #Validating Lambda
        if(is.null(lambda)){
          #not provided by user
          stop("Skewness parameter 'lambda' must be provided for the ESN/SN case.")
        }else
        {
          #validate input
          if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
          if(all(lambda==0)){
            warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
            out = KmomentN(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma)
          }
        }
        if(dist=="SN"){
          out = KmomentESN(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0)
        }else
        {
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else
          {
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = KmomentESN(k = kappa,a = lower,b = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
          }
        }
      }else
      {
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  # cat('\n')
  # call <- match.call()
  # cat("Call:\n")
  # print(call)
  # cat('\n')
  # print(out)
  # cat('\n')
  return(out)
}
# #TESTING
# a = c(-0.8,-0.7,-0.6,-0.5)
# b = c(0.5,0.6,0.7,0.8)
# mu = c(0.1,0.2,0.3,0.4)
# S = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),nrow = length(mu),ncol = length(mu),byrow = TRUE)
#
#
# a = c(-0.8,-0.7,-0.6)
# b = c(0.5,0.6,0.7)
# mu = c(0.1,0.2,0.3)
# Sigma = S = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),nrow = length(mu),ncol = length(mu),byrow = TRUE)
#
# #
# resN = momentsTMD(kappa = c(2,2,2),lower = a,upper = b,mu,varcov = S,dist = "normal")
# #
# resT = momentsTMD(kappa = c(2,2,2),lower = a,upper = b,mu,varcov = S,dist = "t",nu = 8)


# res1 = momentsTMD(kappa = c(2,2),lower = c(0,0),upper = c(Inf,Inf),mus,Ss,dist = "normal",nu=10)
# k = c(5,0)
# mu = mus
# Sigma = Ss
#
# 5
# KmomentT(k = c(2,2),a = c(0,1),b = c(2,4),mu = mus,Sigma = Ss,nu = 200)
# KmomentN(k = c(2,2),a = c(0,1),b = c(2,4),mu = mus,S = Ss)
# momentsTMD(kappa = c(2,2),lower = c(0,0),upper = c(Inf,Inf),mus,Ss)
# KmomentN(k = c(1,1,1),a,b,mu,S)
