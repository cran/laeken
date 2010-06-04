\name{variance}
\alias{variance}
\title{Variance and confidence intervals of Laeken indicators}
\description{
  Compute variance and confidence interval estimates of Laeken indicators.
}
\usage{
variance(inc, weights = NULL, years = NULL, breakdown = NULL, 
         design = NULL, data = NULL, indicator, alpha = 0.05, 
         na.rm = FALSE, type = "bootstrap", \dots)
}
\arguments{
  \item{inc}{either a numeric vector giving the equivalized disposable income, 
    or (if \code{data} is not \code{NULL}) a character string, an integer or a 
    logical vector specifying the corresponding column of \code{data}.}
  \item{weights}{optional; either a numeric vector giving the personal sample 
    weights, or (if \code{data} is not \code{NULL}) a character string, an 
    integer or a logical vector specifying the corresponding column of 
    \code{data}.}
  \item{years}{optional; either a numeric vector giving the different years of 
    the survey, or (if \code{data} is not \code{NULL}) a character string, 
    an integer or a logical vector specifying the corresponding column of 
    \code{data}.  If supplied, values are computed for each year.}
  \item{breakdown}{optional; either a numeric vector giving different strata, 
    or (if \code{data} is not \code{NULL}) a character string, an integer or a 
    logical vector specifying the corresponding column of \code{data}.  If 
    supplied, the values for each stratum are computed in addition to the 
    overall value.}
  \item{design}{optional; either an integer vector or factor giving different 
    strata for stratified sampling designs, or (if \code{data} is not 
    \code{NULL}) a character string, an integer or a logical vector specifying 
    the corresponding column of \code{data}.}
  \item{data}{an optional \code{data.frame}.}
  \item{indicator}{an object inheriting from the class \code{"indicator"} that 
    contains the point estimates of the indicator (see \code{\link{arpr}}, 
    \code{\link{qsr}}, \code{\link{rmpg}} or \code{\link{gini}}).}
  \item{alpha}{a numeric value giving the significance level to be used for 
    computing the confidence interval(s) (i.e., the confidence level is 
    \eqn{1 - }\code{alpha}), or \code{NULL}.}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
  \item{type}{a character string specifying the type of variance estimation 
    to be used.  Currently, only \code{"bootstrap"} is implemented for 
    variance estimation based on bootstrap resampling (see 
    \code{\link{bootVar}}).}
  \item{\dots}{additional arguments to be passed to \code{\link{bootVar}}.}
}
\details{
  This is a wrapper function for computing variance and confidence interval 
  estimates of Laeken indicators.
}
\value{
  An object of the same class as \code{indicator} is returned.  See 
  \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}} or 
  \code{\link{gini}} for details on the components.
}
\author{Andreas Alfons}
\seealso{
\code{\link{bootVar}}, \code{\link{arpr}}, \code{\link{qsr}}, 
\code{\link{rmpg}}, \code{\link{gini}}
}
\examples{
data(eusilc)
a <- arpr("eqIncome", weights = "rb050", data = eusilc)

## naive bootstrap
variance("eqIncome", weights = "rb050", design = "db040", 
    data = eusilc, indicator = a, bootType = "naive", seed = 123)

## bootstrap with calibration
# define auxiliary 0-1 variables for regions
aux <- sapply(levels(eusilc$db040), 
    function(l, x) as.numeric(x == l), 
    x = eusilc$db040)
# retrieve population totals from underlying sample
totals <- sapply(levels(eusilc$db040), 
    function(l, x, w) sum(w[x == l]), 
    x = eusilc$db040, w = eusilc$rb050)
# bootstrap variance
variance("eqIncome", weights = "rb050", design = "db040", 
    data = eusilc, indicator = a, X = aux, totals = totals, 
    seed = 123)
}
\keyword{survey}
