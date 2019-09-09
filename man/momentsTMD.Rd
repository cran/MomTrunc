\name{momentsTMD}
\alias{momentsTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moments for doubly truncated multivariate distributions
}
\description{
It computes the kappa-th raw moment for the doubly truncated \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Student's t-distribution. It also output some other lower moments (than \code{kappa}) involved in the recurrence approach.
}
\usage{
momentsTMD(kappa,lower = NULL,upper = NULL,mu,Sigma,lambda = NULL,tau = NULL,
dist,nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
moments vector of length \eqn{p}. All its elements must be integers greater or equal to \eqn{0}.
}
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
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}.
}
\value{
A data frame containing \eqn{p+2} columns. The \eqn{p} first containing the set of moments involved in the recursive approach and the last two columns containing the \eqn{F} function value (see Galarza and Lachos, 2018) and the expected value. Normal cases (ESN, SN and normal) return \code{prod(kappa)+1} moments while the Student's t-distribution case returns \code{sum(kappa)}+\code{1}. See example section.

Henceforth, we HIGHLY recomend to check the pdf manual instead because of formulae.

The \eqn{F} function is simply

\deqn{F_{\kappa}(\mathbf{a,b},\mu,\Sigma,\nu) = \int_{\mathbf{a}}^{\mathbf{b}} \mathbf{x}^\kappa f(\mathbf{x}) \mathrm{d}\mathbf{x},}
where \eqn{\mathbf{a}} and \eqn{\mathbf{b}} are vectors of length \eqn{p} representing the lower and upper bounds. We have used the short notation \eqn{\mathbf{x}^\kappa = x_1^{\kappa_1} x_2^{\kappa_2}\ldots x_p^{\kappa_p}}. It is easy to see that \eqn{P(\mathbf{a}\le \mathbf{X}\le \mathbf{b})=F_{\mathbf{0}}(\mathbf{a,b},\mu,\Sigma,\nu)}, i.e., the normalizing constant for the doubly truncated density. Then the expected value will be given by \eqn{E[\mathbf{x}^\kappa] = F_\kappa(\mathbf{a,b},\mu,\Sigma,\nu)/F_{\mathbf{0}}(\mathbf{a,b},\mu,\Sigma,\nu)}.

Normal case returns \code{prod(kappa)+1} moments while the Student's t-distribution case returns \code{sum(kappa)}+\code{1}. See example section. For the extended skew-normal case, we recommend to the readers to check the reference.
}

\references{
Kan, R., & Robotti, C. (2017). On Moments of Folded and Truncated Multivariate Normal Distributions. Journal of Computational and Graphical Statistics, (just-accepted).

C.E. Galarza, L.A. Matos, D.K. Dey & V.H. Lachos. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degree of freedom must be a positive integer. If \code{nu >= 200}, Normal case is considered."
}

\section{Warning }{For the Student-t case, the \code{kappa}-\eqn{th} moment can only be computed when \code{sum(kappa)} \eqn{\le} \code{nu}-2.}

\seealso{
\code{\link{meanvarTMD}}, \code{\link{meanvarFMD}},\code{\link{momentsFMD}},\code{\link{meanvarFMD}},\code{\link{dmvESN}},\code{\link{rmvESN}}
}

\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
value1 = momentsTMD(c(2,0,1),a,b,mu,Sigma,dist="normal")
value2 = momentsTMD(c(2,0,1),a,b,mu,Sigma,dist = "t",nu = 7)
value3 = momentsTMD(c(2,0,1),a,b,mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = momentsTMD(c(2,0,1),a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
