\name{drmvSN}
\alias{dmvSN}
\alias{rmvSN}
\title{Multivariate Skew Normal Density and Random Deviates}
\description{
  These functions provide the density function and a random number
  generator for the multivariate skew normal (SN) distribution with mean vector \eqn{mu}, scale matrix \eqn{sigma} and skewness parameter \eqn{lambda}.
}
\usage{
dmvSN(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda)
rmvSN(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda)
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
  a numeric vector of length \eqn{p} representing the skewness parameter for SN and SN cases. If \code{lambda == 0}, the SN/SN reduces to a normal (symmetric) distribution.
  }
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
\code{\link{pmvESN}}, \code{\link{meanvarFMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
}
\examples{
#Univariate case
dmvSN(x = -1,mu = 2,Sigma = 5,lambda = -2)
rmvSN(n = 100,mu = 2,Sigma = 5,lambda = -2)
#Multivariate case
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
lambda = c(-2,0,1,2)
#One observation
dmvSN(x = c(-2,-1,0,1),mu,Sigma,lambda)
rmvSN(n = 100,mu,Sigma,lambda)
#Many observations as matrix
x = matrix(rnorm(4*10),ncol = 4,byrow = TRUE)
dmvSN(x = x,mu,Sigma,lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
