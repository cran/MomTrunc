#MEAN AND VARIANCE

meanvarFMD = function(mu,Sigma,dist = "normal",nu = NULL)
{
  #Validating dims data set
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
    if(is.null(nu)){
      stop("Degrees of freedom 'nu' must be provided for the T case.")}else{
        if(nu < 2){stop("The first moment exists only when the degree of freedom is larger than 2.")
        }else{
          if(nu >= 100){
            warning("For degrees of freedom >= 100, Normal case is considered.",immediate. = TRUE)
            out = meanvarFN(mu = mu,Sigma = Sigma)
          }else{
            if(nu < 4){
              warning("The  second moment exists only when the degree of freedom is larger than 3.",immediate. = TRUE)
              out = meanvarFT(mu = mu,Sigma = Sigma,nu = nu)
            }else{
              out = meanvarFT(mu = mu,Sigma = Sigma,nu = nu)
            }
          }
        }
      }
  }else{
    if(dist != "normal"){stop("The dist values are 'normal' and 't'.")}else{
      if(!is.null(nu)){warning("Nu parameter not considered for normal case.",immediate. = TRUE)}
      nu = 100
      out = meanvarFN(mu = mu,Sigma = Sigma)
    }
  }
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  cat("Mean:\n")
  print(out$mean)
  if(dist == "normal" | (dist == "t" & nu >= 4)){
    cat('\n')
    if(length(mu)==1){cat("Variance:\n")}else{cat("Varcov matrix:\n")}
    print(out$varcov)
    cat('\n')
  }
  return(out)
}

#TESTING
# a = c(-0.8,-0.7,-0.6,-0.5)
# b = c(0.5,0.6,0.7,0.8)
# mu = c(0.1,0.2,0.3,0.4)
# S = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),nrow = length(mu),ncol = length(mu),byrow = TRUE)

# resN = meanvarFMD(mu,S,dist = "normal")
# resT = meanvarFMD(mu,S,dist = "t",nu=4)
