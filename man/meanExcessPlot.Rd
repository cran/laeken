\name{meanExcessPlot}
\alias{meanExcessPlot}
\title{Mean excess plot}
\description{
The Mean Excess plot is a graphical method for detecting the threshold 
(scale parameter) of a Pareto distribution.
}
\usage{
meanExcessPlot(x, w = NULL, probs, interactive = TRUE, ...)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{w}{an optional numeric vector giving sample weights.}
  \item{probs}{an optional numeric vector of probabilities with values in 
    \eqn{[0,1]}, defining the quantiles to be plotted.  This is useful for 
    large data sets, when it may not be desirable to plot every single point.}
  \item{interactive}{a logical indicating whether the threshold (scale 
    parameter) can be selected interactively by clicking on points.  
    Information on the selected threshold is then printed on the \R console.}
  \item{\dots}{additional arguments to be passed to 
    \code{\link[graphics]{plot.default}}.}
}
\details{
The corresponding mean excesses are plotted against the values of \code{x} (if 
supplied, only those specified by \code{probs}).  If the tail of the data 
follows a Pareto distribution, these observations show a positive linear trend.
The leftmost point of a fitted line can thus be used as an estimate of the 
threshold (scale parameter).

The interactive selection of the threshold (scale parameter) is implemented 
using \code{\link[graphics]{identify}}.  For the usual \code{X11} device, the 
selection process is thus terminated by pressing any mouse button other than 
the first.  For the \code{quartz} device (on Mac OS X systems), the process is 
terminated either by a secondary click (usually second mouse button or 
\code{Ctrl}-click) or by pressing the \code{ESC} key.
}
\value{
If \code{interactive} is \code{TRUE}, the last selection for the threshold is 
returned invisibly as an object of class \code{"paretoScale"}, which consists 
of the following components:
\item{x0}{the selected threshold (scale parameter).}
\item{k}{the number of observations in the tail (i.e., larger than the 
  threshold).}
}
\note{
  The functionality to account for sample weights and to select the threshold 
  (scale parameter) interactively was introduced in version 0.2.  
}
\author{Andreas Alfons and Josef Holzer}
\seealso{
  \code{\link{minAMSE}}, \code{\link{paretoScale}}, \code{\link{paretoQPlot}}, 
  \code{\link[graphics]{identify}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
eusilc <- eusilc[!duplicated(eusilc$db030),]

# with sample weights
meanExcessPlot(eusilc$eqIncome, w = eusilc$db090)

# without sample weights
meanExcessPlot(eusilc$eqIncome)
}
\keyword{hplot}
