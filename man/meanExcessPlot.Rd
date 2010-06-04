\name{meanExcessPlot}
\alias{meanExcessPlot}
\title{Mean excess plot}
\description{
The Pareto quantile plot is a graphical method for detecting the threshold 
(scale parameter) of a Pareto distribution.
}
\usage{
meanExcessPlot(x, probs, ...)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{probs}{a numeric vector of probabilities with values in \eqn{[0,1]}, 
    defining the quantiles to be plotted.  This is useful for large data sets, 
    when it may not be desirable to plot every single point.}
  \item{\dots}{additional arguments to be passed to 
    \code{\link[graphics]{plot.default}}.}
}
\details{
The corresponding mean excesses are plotted against the values of \code{x} (if 
supplied, only those specified by \code{probs}).  If the tail of the data 
follows a Pareto distribution, these observations show a positive linear trend.
The leftmost point of a fitted line can thus be used as an estimate of the 
scale parameter.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{paretoQPlot}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
meanExcessPlot(eusilc$eqIncome[!duplicated(eusilc$db030)])
}
\keyword{hplot}
