\name{meanvarTMD}
\alias{meanvarTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and variance for doubly truncated multivariate distributions
}
\description{
It computes the mean vector and variance-covariance matrix for the doubly truncated \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Student's t-distribution.
}
\usage{
meanvarTMD(lower = NULL,upper = NULL,mu,Sigma,lambda = NULL,tau = NULL,dist,nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
a numeric vector of length \eqn{p} representing the skewness parameter for SN and ESN cases. If \code{lambda == 0}, the ESN/SN reduces to a normal (symmetric) distribution.
}
  \item{tau}{
It represents the extension parameter for the ESN distribution. If \code{tau == 0}, the ESN reduces to a SN distribution.
}
  \item{dist}{
represents the folded distribution to be computed. The values are \code{normal}, \code{SN} , \code{ESN} and \code{t}  for the doubly truncated Normal, Skew-normal, Extended Skew-normal and Student's t-distribution respectively.
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution.
}
}
\details{
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}. Normal case code is an R adaptation of the Matlab available function \code{dtmvnmom.m} from Kan & Robotti (2017) and it is used for \code{p<=3}. For higher dimensions we use an extension of the algorithm in Vaida (2009).
}
\value{
It returns a list with three elements:
\item{mean}{the mean vector of length \eqn{p}}
\item{EYY}{the second moment matrix of dimensions \eqn{p}x\eqn{p}}
\item{varcov}{the variance-covariance matrix of dimensions \eqn{p}x\eqn{p}}
}
\references{
Kan R. & Robotti C. (2017) On Moments of Folded and Truncated Multivariate Normal Distributions, Journal of Computational and Graphical Statistics, 26:4, 930-934.

C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.

Vaida, F. & Liu, L. (2009). Fast implementation for normal mixed effects models with censored response. Journal of Computational and Graphical Statistics, 18(4), 797-817.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degree of freedom must be a positive integer. If \code{nu >= 200}, Normal case is considered."
}

\section{Warning }{The mean can only be provided when \code{nu} is larger than 2. By the other hand, the varcov matrix can only be provided when \code{nu} is larger than 3.}

\seealso{
\code{\link{momentsTMD}}, \code{\link{meanvarFMD}}, \code{\link{momentsFMD}},\code{\link{dmvESN}},\code{\link{rmvESN}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
value1 = meanvarTMD(a,b,mu,Sigma,dist="normal")
value2 = meanvarTMD(a,b,mu,Sigma,dist = "t",nu = 4)
value3 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
