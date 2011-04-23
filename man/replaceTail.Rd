\name{replaceTail}
\alias{replaceTail}
\alias{replaceTail.paretoTail}
\alias{replaceOut}
\title{Replace observations under a Pareto model}
\description{
  Replace observations under a Pareto model for the upper tail with values 
  drawn from the fitted distribution.
}
\usage{
replaceTail(x, \dots)

\method{replaceTail}{paretoTail}(x, all = TRUE, \dots)

replaceOut(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"paretoTail"} (see 
    \code{\link{paretoTail}}).}
  \item{all}{a logical indicating whether all observations in the upper tail 
    should be replaced or only those flagged as outliers.}
  \item{\dots}{additional arguments to be passed down.}
}
\details{
  \code{replaceOut(x, \dots)} is a simple wrapper for \code{replaceTail(x, 
  all = FALSE, \dots)}.
}
\value{
  A numeric vector consisting mostly of the original values, but with 
  observations in the upper tail replaced with values from the fitted Pareto 
  distribution.
}
\author{Andreas Alfons}
\seealso{
  \code{\link{paretoTail}}, \code{\link{reweightOut}}, \code{\link{shrinkOut}} 
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

# replacement of outliers
eqIncome <- replaceOut(fit)
gini(eqIncome, weights = eusilc$rb050)

# replacement of whole tail
eqIncome <- replaceTail(fit)
gini(eqIncome, weights = eusilc$rb050)
}
\keyword{manip}
