\name{paretoScale}
\alias{paretoScale}
\title{Estimate the scale parameter of a Pareto distribution}
\description{
  Estimate the scale parameter of a Pareto distribution, i.e., the threshold 
  for Pareto tail modeling.
}
\usage{
paretoScale(x, w = NULL, groups = NULL, method = "vanKerm", na.rm = FALSE)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{w}{an optional numeric vector giving sample weights.}
  \item{groups}{an optional vector or factor specifying groups of elements of 
    \code{x} (e.g., households).  If supplied, each group of observations is 
    expected to have the same value in \code{x} (e.g., household income).  Only 
    the values of every first group member to appear are used for estimating 
    the threshold (scale parameter).}
  \item{method}{a character string specifying the estimation method.  If 
    \code{"vanKerm"}, van Kerm's method is used, which is a rule of thumb 
    specifically designed for the equivalized disposable income in EU-SILC data 
    (currently the only method implemented).}
  \item{na.rm}{a logical indicating whether missing values in \code{x} should 
    be omitted.}
}
\value{
An object of class \code{"paretoScale"}, which consists of the following 
components:
\item{x0}{the threshold (scale parameter).}
\item{k}{the number of observations in the tail (i.e., larger than the 
  threshold).}
}
\author{Andreas Alfons}
\references{
Van Kerm, P. (2007) Extreme incomes and the estimation of poverty and 
inequality indicators from EU-SILC. IRISS Working Paper Series 2007-01, 
CEPS/INSTEAD.
}
\seealso{
  \code{\link{minAMSE}}, \code{\link{paretoQPlot}}, \code{\link{meanExcessPlot}}
}
\examples{
data(eusilc)
paretoScale(eusilc$eqIncome, eusilc$db090, groups = eusilc$db030)
}
\keyword{manip}
