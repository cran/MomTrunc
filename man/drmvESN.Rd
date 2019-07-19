\name{drmvESN}
\alias{dmvESN}
\alias{rmvESN}
\title{Multivariate Extended-Skew Normal Density and Random Deviates}
\description{
  These functions provide the density function and a random number
  generator for the multivariate extended-skew normal (ESN) and skew-normal (SN, as particular case) distribution with mean vector \eqn{mu}, scale matrix \eqn{sigma}, skewness parameter \eqn{lambda} and extension parameter \eqn{tau}.
}
\usage{
dmvESN(x,mu=c(0,0),Sigma=diag(2),lambda=c(-1,1),tau=1)
rmvESN(n,mu=c(0,0),Sigma=diag(2),lambda=c(-1,1),tau=1)
}
\arguments{
   \item{x}{vector or matrix of quantiles. If \code{x} is a matrix, each
     row is taken to be a quantile.}
   \item{n}{number of observations.}
   \item{mu}{
  a numeric vector of length \eqn{p} representing the location parameter.
  }
  \item{Sigma}{
  a numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the scale parameter.
  }
  \item{lambda}{
  a numeric vector of length \eqn{p} representing the skewness parameter for SN and ESN cases. If \code{lambda == 0}, the ESN/SN reduces to a normal (symmetric) distribution.
  }
  \item{tau}{
  It represents the extension parameter for the ESN distribution. If \code{tau == 0}, the ESN reduces to a SN distribution.
  }
}

\details{
The \code{rmvESN} function is based on a rejection algorithm, so for small values of \code{tau}, that is, when \eqn{\tau \le -2\sqrt{1+\sum{\lambda_i^2}}}, it is possible that the algorithm return less than the \code{n} observations and also take more time than usual, due to the high rejection rate.
}
\references{
C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pmvESN}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
}
\examples{
#Univariate case
dmvESN(x = -1,mu = 2,Sigma = 5,lambda = -2,tau = 0.5)
rmvESN(n = 100,mu = 2,Sigma = 5,lambda = -2,tau = 0.5)
#Multivariate case
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
lambda = c(-2,0,1,2)
tau = 2
#One observation
dmvESN(x = c(-2,-1,0,1),mu,Sigma,lambda,tau)
rmvESN(n = 100,mu,Sigma,lambda,tau)
#Many observations as matrix
x = matrix(rnorm(4*10),ncol = 4,byrow = TRUE)
dmvESN(x = x,mu,Sigma,lambda,tau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
