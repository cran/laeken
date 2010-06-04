\name{thetaTM}
\alias{thetaTM}
\title{Trimmed mean estimator}
\description{
  Estimate the shape parameter of a Pareto distribution using a trimmed mean 
  approach.
}
\usage{
thetaTM(x, k, beta = 0.05)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{the number of observations in the upper tail to which the Pareto 
    distribution is fitted.}
  \item{beta}{A numeric vector of length two giving the trimming proportions 
    for the lower and upper end of the tail, respectively.  If a single numeric 
    value is supplied, it is recycled.}
}
\value{
  The estimated shape parameter.
}
\references{
Brazauskas, V. and Serfling, R. (2000) Robust estimation of tail parameters 
for two-parameter Pareto and exponential models via generalized quantile 
statistics. \emph{Extremes}, \bold{3}(3), 231--249.

Brazauskas, V. and Serfling, R. (2000) Robust and efficient estimation of the 
tail index of a single-parameter Pareto distribution. \emph{North American 
Actuarial Journal}, \bold{4}(4), 12--27.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaTM(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
