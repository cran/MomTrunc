#FULL MOMENTS FOLDED

momentsFMD = function(kappa,mu,Sigma,dist = "normal",nu = NULL)
{
  #Validating dims data set
  if(ncol(as.matrix(kappa)) > 1 | !all(kappa >= 0) | length(c(kappa)) != length(c(mu))) stop("kappa must be numeric with same dimensions than mu.")
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("mu must be numeric and have just one column")

  #validate mean an Sigma dimensions

  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.positive.definite(Sigma))stop("Sigma must be a square symmetrical real posite definite matrix.")
  }

  #validating distributions and nu parameter
  if(dist == "t"){
    if(is.null(nu) | nu%%1!=0){
      stop("Degrees of freedom 'nu' must be a positive integer provided for the T case.")}else{
        if(!all(kappa == 0) & nu < max(3,sum(kappa)+2)){stop("The kappa-th moment exists only when the degree of freedom is greater than or equal to 'sum(kappa)+2'.")}
        if(nu >= 100){
          warning("For degrees of freedom >= 100, Normal case is considered.",immediate. = TRUE)
          out = KmomentFN(k = kappa,mu = mu,Sigma = Sigma)
        }else{
          out = KmomentFT(k = kappa,mu = mu,Sigma = Sigma,nu = nu)
        }
      }
  }else{
    if(dist != "normal"){stop("The dist values are 'normal' and 't'.")}else{
      if(!is.null(nu)){warning("Nu parameter not considered for normal case.",immediate. = TRUE)}
      out = KmomentFN(k = kappa,mu = mu,Sigma = Sigma)
    }
  }
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  print(out)
  cat('\n')
  return(out)
}
