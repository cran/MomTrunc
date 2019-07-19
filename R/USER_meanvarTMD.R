#MEAN AND VARIANCE

meanvarTMD = function(lower = NULL,upper = NULL,mu,Sigma,lambda = NULL, tau = NULL,dist,nu = NULL)
{
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
    out = meanvarN7(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
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
              out = meanvarN7(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
            }else{
              if(nu < 4){
                warning("Sorry, we can only compute the second moment when the degrees of freedom is larger than 3.",immediate. = TRUE)
                out = meanvarT(a = lower,b = upper,mu = mu,Sigma = Sigma,nu = nu)
              }else{
                out = meanvarT(a = lower,b = upper,mu = mu,Sigma = Sigma,nu = nu)
              }
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
            out = meanvarN7(lower = lower,upper = upper,mu = mu,Sigma = Sigma)
          }
        }
        if(dist=="SN"){
          out = meanvarESN7(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0)
        }else{
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else{
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = meanvarESN7(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
          }
        }
      }else{
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  return(out)
}

# #TESTING
# a = c(-0.8,-0.7,-0.6,-0.5)
# b = c(0.5,0.6,0.7,0.8)
# mu = c(0.1,0.2,0.3,0.4)
# S = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),nrow = length(mu),ncol = length(mu),byrow = TRUE)
#
#
# #NORMAL 0.185
# ptm <- proc.time()
# for(i in 1:10){
#   resN = meanvarTMD(lower = a,upper = b,mu,S)}
# (proc.time() - ptm)/10
#
# #t RECURRENCE 1.60
# ptm <- proc.time()
# for(i in 1:10){
# resT = meanvarTMD(lower = a,upper = b,mu,S,dist = "t",nu = 5)}
# (proc.time() - ptm)/10
#
# #t SLICE SAMPLING 5.063
# ptm <- proc.time()
# for(i in 1:10){
# resT2 = TT.moment(mu = mu,S = S,nu = 5,lower = a,upper = b)}
# (proc.time() - ptm)/10
#
# library(tmvtnorm)
# ptm <- proc.time()
# gent = rtmvt(n = 10^5,mean = mu,sigma = S,df = 5,lower = a,upper = b,algorithm = "gibbs")
# proc.time() - ptm
# mapprox = apply(gent,2,mean)
# varapprox = var(gent)
#
# round(cbind(resT$mu,resT2$EX,mapprox),4)
# round(cbind(diag(resT$varcov),diag(resT2$CovX),diag(varapprox)),4)
#
#
#
# #p = 10
#
# p = 3 #nosso 88.30
# a  = seq(-0.9,0,length.out = p)
# b  = seq(0,0.9,length.out = p) + seq(0.1,0.6,length.out = p)
# mu = seq(0.25,0.75,length.out = p)*b
# s  = matrix(0.5*rnorm(p^2),p,p)
# Sigma = S = s%*%t(s)
# #S[1,2] = S[2,1] = sqrt((S[1,1]*S[2,2]))
# nu = 5
# #
# # #t RECURRENCE 1.60
# # ptm <- proc.time()
# resT = meanvarTMD(lower = a,upper = b,mu,S)

# proc.time() - ptm
#
# #t SLICE SAMPLING 5.063
# ptm <- proc.time()
# resT2 = TT.moment(mu = mu,S = S,nu = nu,lower = a,upper = b)
# proc.time() - ptm
#
# round(cbind(resT$mu,resT2$EX),4)
# round(cbind(diag(resT$varcov),diag(resT2$CovX)),4)
#
#
# #nosso =
# #lin =
#
#
#
#
#
#
# library(tmvtnorm)
# ptm <- proc.time()
# gent = rtmvt(n = 10^5,mean = mu,sigma = S,df = nu,lower = a,upper = b,algorithm = "gibbs")
# proc.time() - ptm
# mapprox = apply(gent,2,mean)
# varapprox = var(gent)
#
# round(cbind(resT$mu,resT2$EX,mapprox),4)
# round(cbind(diag(resT$varcov),diag(resT2$CovX),diag(varapprox)),4)
#
# aux3 = resT$mu
#
# #
# #
# #
# # GB = GenzBretz(maxpts = 5e4, abseps = 1e-6, releps = 0)#9.83
# # GB = GenzBretz(maxpts = 5e4, abseps = 1e-5, releps = 0)#9.77
# # GB = GenzBretz(maxpts = 5e4, abseps = 1e-4, releps = 0)#2.00
# # ptm <- proc.time()
# # for(i in 1:20){
# # F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = S, algorithm = GB)[1]}
# # proc.time() - ptm
# #
# #
# #
# # GB = GenzBretz(maxpts = 4e4, abseps = 1e-6, releps = 0)#9.75
# # GB = GenzBretz(maxpts = 4e4, abseps = 1e-5, releps = 0)#9.78
# # GB = GenzBretz(maxpts = 5e4, abseps = 5e-5, releps = 0)#4s
# # ptm <- proc.time()
# # for(i in 1:20){
# #   F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = S, algorithm = GB)[1]}
# # proc.time() - ptm