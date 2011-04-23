\name{shrinkOut}
\alias{shrinkOut}
\alias{shrinkOut.paretoTail}
\title{Shrink outliers in the Pareto model}
\description{
  Shrink observations that are flagged as outliers in a Pareto model for the 
  upper tail of the distribution to the theoretical quantile used for outlier 
  detection.
}
\usage{
shrinkOut(x, \dots)

\method{shrinkOut}{paretoTail}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"paretoTail"} (see 
    \code{\link{paretoTail}}).}
  \item{\dots}{additional arguments to be passed down (currently ignored as 
    there are no additional arguments in the only method implemented).}
}
\value{
  A numeric vector consisting mostly of the original values, but with 
  outlying observations in the upper tail shrunken to the corresponding 
  theoretical quantile of the fitted Pareto distribution.
}
\author{Andreas Alfons}
\seealso{
  \code{\link{paretoTail}}, \code{\link{reweightOut}}, 
  \code{\link{replaceOut}}, \code{\link{replaceTail}}
}
\examples{
data(eusilc)

## gini coefficient without Pareto tail modeling
gini("eqIncome", weights = "rb050", data = eusilc)

## gini coefficient with Pareto tail modeling
# estimate threshold
ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090, 
    groups = eusilc$db030)
# estimate shape parameter
fit <- paretoTail(eusilc$eqIncome, k = ts$k, 
    w = eusilc$db090, groups = eusilc$db030)
# shrink outliers
eqIncome <- shrinkOut(fit)
gini(eqIncome, weights = eusilc$rb050)
}
\keyword{manip}
