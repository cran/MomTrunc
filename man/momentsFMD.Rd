\name{momentsFMD}
\alias{momentsFMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moments for folded multivariate distributions
}
\description{
It computes the kappa-th raw moment for the folded \code{p}-variate Normal and Student's t-distribution. It also output some other lower moments (than \code{kappa}) involved in the recurrence approach.
}
\usage{
momentsFMD(kappa,mu,Sigma,dist = "normal",nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
moments vector of length \eqn{p}. All its elements must be integers greater or equal to \eqn{0}.
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
A data frame containing \eqn{p+1} columns. The \eqn{p} first containing the set of moments involved in the recursive approach and the last column containing the expected value.

Normal case returns \code{prod(kappa)+1} moments while the Student's t-distribution case returns \code{sum(kappa)+1}. See example section.
}

\references{
Kan, R., & Robotti, C. (2017). On Moments of Folded and Truncated Multivariate Normal Distributions. Journal of Computational and Graphical Statistics, (just-accepted).

Galarza C.E. & Lachos, V.H. (2018). On moments of folded and truncated multivariate Student-t distribution: A recurrence approach.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degree of freedom must be a positive integer. If \code{nu >= 100}, Normal case is considered."
}

\section{Warning }{The \code{kappa}-\eqn{th} moment can only be computed when \code{sum(kappa)} \eqn{\le} \code{nu}-2.}

\seealso{
\code{\link{meanvarFMD}},\code{\link{momentsTMD}},\code{\link{meanvarTMD}},\code{\link{cdfFMD}}}

\examples{
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
nu = 7
value1 = momentsFMD(c(2,0,1),mu,Sigma) #normal case
value2 = momentsFMD(c(0,2,0),mu,Sigma,dist = "t",nu) #t case
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Folded}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}