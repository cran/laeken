\name{thetaMoment}
\alias{thetaMoment}
\title{Moment estimator}
\description{
  Estimate the shape parameter of a Pareto distribution based on moments.
}
\usage{
thetaMoment(x, k)
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
Dekkers, A.L.M., Einmahl, J.H.J. and de Haan, L. (1989) A moment estimator for 
the index of an extreme-value distribution. \emph{The Annals of Statistics}, 
\bold{17}(4), 1833--1855.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaMoment(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
