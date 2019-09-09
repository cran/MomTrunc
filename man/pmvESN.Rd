\name{pmvESN}
\alias{pmvESN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Extended-Skew Normal Probabilities
}
\description{
Computes the distribution function of the multivariate extended-skew normal (ESN) and skew-normal (SN, as particular case) distribution for arbitrary limits, mean vector \eqn{mu}, scale matrix \eqn{sigma}, skewness parameter \eqn{lambda} and extension parameter \eqn{tau}.
}
\usage{
pmvESN(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),
        mu = rep(0,length(lambda)),Sigma,lambda,tau, ...)
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
\item{...}{
Arguments passed to the \code{pmvnorm} fucntion from the \code{mvtnorm}, for instance, \code{algorithm}.
}
}
\details{
This code uses the \code{pmvnorm} function from \code{mvtnorm} package. For the normal case (all lambda and tau equal to zero), it is better to use \code{pmvnorm} directly.
}
\value{
The evaluated distribution function is returned.
}
\references{
Genz, A. and Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195. Springer-Verlag, Heidelberg.

C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{dmvESN}},\code{\link{rmvESN}}, \code{\link{meanvarFMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
}
\examples{
lower = rep(-Inf,4)
upper = c(-1,0,2,5)
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
nrow = length(mu),ncol = length(mu),byrow = TRUE)
lambda = c(-2,0,1,2)
tau = 2
pmvESN(lower,upper,mu,Sigma,lambda,tau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
