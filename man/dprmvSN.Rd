\name{dprmvSN}
\alias{dmvSN}
\alias{pmvSN}
\alias{rmvSN}
\title{Multivariate Skew Normal Density and Probabilities and Random Deviates}
\description{
  These functions provide the density function and a random number
  generator for the multivariate skew normal (SN) distribution with mean vector \code{mu}, scale matrix \code{Sigma} and skewness parameter \code{lambda}.
}
\usage{
dmvSN(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda)
pmvSN(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),
        mu = rep(0,length(lambda)),Sigma,lambda,log2 = FALSE)
rmvSN(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda)
}
\arguments{
   \item{x}{vector or matrix of quantiles. If \code{x} is a matrix, each
     row is taken to be a quantile.}
   \item{n}{number of observations.}
    \item{lower}{
the vector of lower limits of length \eqn{p}.
}
\item{upper}{
the vector of upper limits of length \eqn{p}.
}
   \item{mu}{
  a numeric vector of length \eqn{p} representing the location parameter.
  }
  \item{Sigma}{
  a numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the scale parameter.
  }
  \item{lambda}{
  a numeric vector of length \eqn{p} representing the skewness parameter for SN and SN cases. If \code{lambda == 0}, the SN/SN reduces to a normal (symmetric) distribution.
  }
      \item{log2}{
  a boolean variable, indicating if the log2 result should be returned. This is useful when the true probability is too small for the machine precision.
  }
}

\value{
\code{dmvSN} gives the density, \code{pmvSN} gives the distribution function, and \code{rmvSN} generates random deviates for the Multivariate Skew-normal Distribution.
}

\references{
Galarza, C. E., Matos, L. A., Dey, D. K., & Lachos, V. H. (2022a). "On moments of folded and doubly truncated multivariate extended skew-normal distributions." Journal of Computational and Graphical Statistics, 1-11 <doi:10.1080/10618600.2021.2000869>.

Galarza, C. E., Matos, L. A., Castro, L. M., & Lachos, V. H. (2022b). Moments of the doubly truncated selection elliptical distributions with emphasis on the unified multivariate skew-t distribution. Journal of Multivariate Analysis, 189, 104944 <doi:10.1016/j.jmva.2021.104944>.

Galarza, C.E., Matos, L.A. and Lachos, V.H. (2022c). An EM algorithm for estimating the parameters of the multivariate skew-normal distribution with censored responses. Metron. <doi:10.1007/s40300-021-00227-4>.

Genz, A.,  (1992) "Numerical computation of multivariate normal probabilities," Journal of Computational and Graphical Statistics, 1, 141-149 <doi:10.1080/10618600.1992.10477010>.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{dmvESN}}, \code{\link{pmvESN}}, \code{\link{rmvESN}}, \code{\link{meanvarFMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
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

lower = rep(-Inf,4)
upper = c(-1,0,2,5)
pmvSN(lower,upper,mu,Sigma,lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
