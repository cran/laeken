\name{thetaHill}
\alias{thetaHill}
\title{Hill estimator}
\description{
  The Hill estimator uses the maximum likelihood principle to estimate the 
  shape parameter of a Pareto distribution.
}
\usage{
thetaHill(x, k)
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
Hill, B.M. (1975) A simple general approach to inference about the tail of a 
distribution. \emph{The Annals of Statistics}, \bold{3}(5), 1163--1174.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}, \code{\link{thetaPDC}}, \code{\link{thetaWML}}, 
  \code{\link{thetaISE}}, \code{\link{minAMSE}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaHill(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
