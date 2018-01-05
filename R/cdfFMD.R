cdfFMD = function(x,mu,Sigma,dist = "normal",nu = NULL)
{
  #Validating dims data set
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("y must be numeric and have just one column")
  if(ncol(as.matrix(x)) > 1 | !all(x >= 0) | length(c(x)) != length(c(mu))) stop("x must be numeric with same dimensions than mu.")

  #validate mean an Sigma dimensions

  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.positive.definite(Sigma))stop("Sigma must be a square symmetrical real posite definite matrix.")
  }

  #validating distributions and nu parameter
  if(dist == "t"){
    if(is.null(nu)){
      stop("Degrees of freedom 'nu' must be provided for the T case.")}else{
        if(nu <= 0){stop("Degree of freedom must be a positive integer.")
        }else{
          if(nu >= 100){
            warning("For degrees of freedom >= 100, Normal case is considered.",immediate. = TRUE)
            out = cdfFT(x = x,mu = mu,Sigma = Sigma,nu = 0)
          }else{
            out = cdfFT(x = x,mu = mu,Sigma = Sigma,nu = nu)
          }
        }
      }
  }else{
    if(dist != "normal"){stop("The dist values are 'normal' and 't'.")}else{
      if(!is.null(nu)){warning("Nu parameter not considered for normal case.",immediate. = TRUE)}
      out = cdfFT(x = x,mu = mu,Sigma = Sigma,nu = 0)
    }
  }
  return(out)
}
