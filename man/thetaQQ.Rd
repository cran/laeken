\name{thetaQQ}
\alias{thetaQQ}
\title{QQ-estimator}
\description{
  Estimate the shape parameter of a Pareto distribution using a 
  quantile-quantile approach.
}
\usage{
thetaQQ(x, k)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{the number of observations in the upper tail to which the Pareto 
    distribution is fitted.}
}
\value{
  The estimated shape parameter.
}
\references{
Kratz, M.F. and Resnick, S.I. (1996) The QQ-estimator and heavy tails. 
\emph{Stochastic Models}, \bold{12}(4), 699--724.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaQQ(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
