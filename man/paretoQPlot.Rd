\name{paretoQPlot}
\alias{paretoQPlot}
\title{Pareto quantile plot}
\description{
The Pareto quantile plot is a graphical method for inspecting the parameters of 
a Pareto distribution.
}
\usage{
paretoQPlot(x, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{\dots}{additional arguments to be passed to 
    \code{\link[graphics]{plot.default}}.}
}
\details{
The logarithms of the observed values are plotted against the corresponding 
theoretical quantiles of the standard exponential distribution.  If the tail 
of the data follows a Pareto distribution, these observations form (almost) a 
straight line.  The leftmost point of a fitted line can thus be used as an 
estimate of the scale parameter. The slope of the fitted line is in turn an 
estimate of \eqn{\frac{1}{\theta}}{1/theta}, the reciprocal of the shape 
parameter.
}
\references{
Beirlant, J., Vynckier, P. and Teugels, J.L. (1996) Tail index estimation, 
Pareto quantile plots, and regression diagnostics. \emph{Journal of the 
American Statistical Association}, \bold{91}(436), 1659--1667.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{minAMSE}}, \code{\link{meanExcessPlot}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
paretoQPlot(eusilc$eqIncome[!duplicated(eusilc$db030)])
}
\keyword{hplot}
