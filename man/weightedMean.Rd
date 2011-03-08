\name{weightedMean}
\alias{weightedMean}
\title{Weighted mean}
\description{
  Compute the weighted mean.
}
\usage{
weightedMean(x, weights = NULL, na.rm = FALSE)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{weights}{an optional numeric vector giving the sample weights.}
  \item{na.rm}{a logical indicating whether missing values in \code{x} should 
    be omitted.}
}
\details{
  This is a simple wrapper function calling \code{\link[stats]{weighted.mean}} 
  if sample weights are supplied and \code{\link{mean}} otherwise.
}
\value{
  The weighted mean of values in \code{x} is returned.
}
\author{Andreas Alfons}
\seealso{\code{\link{incMean}}}
\examples{
data(eusilc)
weightedMean(eusilc$eqIncome, eusilc$rb050)
}
\keyword{survey}
