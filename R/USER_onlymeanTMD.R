#TN_mean
onlymeanTMD = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma,lambda = NULL, tau = NULL,
                       dist, nu = NULL){
  mu = c(mu)
  lambda = c(lambda)
  #Validating dims data set
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
    out = onlymeanN(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
  }else{
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else{
        if(nu%%1!=0){
          stop("Degrees of freedom 'nu' must be an integer greater than 2.")
        }else{
          if(nu <= 2){stop("Sorry, we can only compute the first moment for degrees of freedom larger than 2.")
          }else{
            if(nu >= 200){
              warning("For degrees of freedom >= 200, Normal case is considered.",immediate. = TRUE)
              out = onlymeanN(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
            }else{
              out = onlymeanT(a = lower,b = upper,mu = mu,Sigma = Sigma,nu = nu)
            }
          }
        }
      }
    }else{
      if(dist == "ESN" | dist == "SN"){
        #Validating Lambda
        if(is.null(lambda)){
          #not provided by user
          stop("Skewness parameter 'lambda' must be provided for the ESN/SN case.")
        }else{
          #validate input
          if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
          if(all(lambda==0)){
            warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
            out = onlymeanN(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
          }
        }
        if(dist=="SN"){
          out = onlymeanESN(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0)
        }else{
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else{
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = onlymeanESN(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
          }
        }
      }else{
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  return(as.numeric(unlist(out)))
}
