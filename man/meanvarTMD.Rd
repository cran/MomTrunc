\name{meanvarTMD}
\alias{meanvarTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and variance for doubly truncated multivariate distributions
}
\description{
It computes the mean vector and variance-covariance matrix for the doubly truncated \code{p}-variate Normal and Student's t-distribution.
}
\usage{
meanvarTMD(lower = NULL,upper = NULL,mu,Sigma,dist = "normal",nu = NULL)
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
  \item{dist}{
represents the truncated distribution to be computed. The values are \code{normal} (by default) for the doubly truncated Normal distribution and \code{t} for the truncated Student's t-distribution.
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution.
}
}
\details{
Normal case by default, i.e., when \code{dist} is not provided. Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}. Normal case code is an R adaptation of the Matlab available function \code{dtmvnmom.m} from Kan & Robotti (2017).
}
\value{
It returns a list with two elements:
\item{mean}{the mean vector of length \eqn{p}}
\item{varcov}{the variance-covariance matrix of dimensions \eqn{p}x\eqn{p}}
}
\references{
Kan, R., & Robotti, C. (2017). On Moments of Folded and Truncated Multivariate Normal Distributions. Journal of Computational and Graphical Statistics, (just-accepted).

Galarza C.E. & Lachos, V.H. (2018). On moments of folded and doubly truncated multivariate Student-t distribution: A recurrence approach.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degree of freedom must be a positive integer. If \code{nu >= 100}, Normal case is considered."
}

\section{Warning }{The mean can only be computed when \code{nu} is larger than 2. By the other hand, the varcov matrix can only be computed when \code{nu} is larger than 3.}

\seealso{
\code{\link{momentsTMD}}, \code{\link{meanvarFMD}}, \code{\link{momentsFMD}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
nu = 4
value1 = meanvarTMD(a,b,mu,Sigma) #normal case
value2 = meanvarTMD(a,b,mu,Sigma,dist = "t",nu) #t case
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
