\name{cdfFMD}
\alias{cdfFMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative distribution function for folded multivariate distributions
}
\description{
It computes the cumulative distribution function on \code{x} for a folded \code{p}-variate Normal and Student's t-distribution.
}
\usage{
cdfFMD(x,mu,Sigma,dist = "normal",nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of length \eqn{p} where the cdf is evaluated.
}
  \item{mu}{
a numeric vector of length \eqn{p} representing the location parameter.
}
  \item{Sigma}{
a numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the scale parameter.
}
  \item{dist}{
represents the folded distribution to be computed. The values are \code{normal} (by default) for the folded Normal distribution and \code{t} for the folded Student's t-distribution.
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution.
}
}
\details{
Normal case by default, i.e., when \code{dist} is not provided. Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}.
}
\value{
It returns the distribution value for a single point \code{x}.
}
\references{
Psarakis, S., & Panaretos, J. (2001). On some bivariate extensions of the folded normal and the folded-t distributions.

Chakraborty, A. K., & Chatterjee, M. (2013). On multivariate folded normal distribution. Sankhya B, 75(1), 1-15.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degree of freedom must be a positive integer. If \code{nu >= 100}, Normal case is considered."
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{momentsFMD}}, \code{\link{meanvarFMD}}
}
\examples{
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
nu = 4
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma) #normal case
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma,dist = "t",nu) #t case
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Folded}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Student's t}
