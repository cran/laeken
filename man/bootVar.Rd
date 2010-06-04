\name{bootVar}
\alias{bootVar}
\title{Bootstrap variance and confidence intervals of Laeken indicators}
\description{
  Compute variance and confidence interval estimates of Laeken indicators based 
  on bootstrap resampling.
}
\usage{
bootVar(inc, weights = NULL, years = NULL, breakdown = NULL, 
        design = NULL, data = NULL, indicator, R = 100, 
        bootType = c("calibrate", "naive"), X, totals, 
        ciType = c("perc", "norm", "basic"),
        alpha = 0.05, seed = NULL, na.rm = FALSE, ...)
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
    the corresponding column of \code{data}.  If supplied, this is used as 
    \code{strata} argument in the call to \code{\link[boot]{boot}}.}
  \item{data}{an optional \code{data.frame}.}
  \item{indicator}{an object inheriting from the class \code{"indicator"} that 
    contains the point estimates of the indicator (see \code{\link{arpr}}, 
    \code{\link{qsr}}, \code{\link{rmpg}} or \code{\link{gini}}).}
  \item{R}{a numeric value giving the number of bootstrap replicates.}.
  \item{bootType}{a character string specifying the type of bootstap to be 
    performed.  Possible values are \code{"calibrate"} (for calibration of the 
    sample weights of the resampled observations in every iteration) and 
    \code{"naive"} (for a naive bootstrap without calibration of the sample 
    weights).}
  \item{X}{if \code{bootType} is \code{"calibrate"}, a matrix of calibration 
    variables.}
  \item{totals}{numeric; if \code{bootType} is \code{"calibrate"}, this gives 
    the population totals.  If \code{years} is \code{NULL}, a vector should be 
    supplied, otherwise a matrix in which each row contains the population 
    totals of the respective year.}
  \item{ciType}{a character string specifying the type of confidence interval(s) 
    to be computed.  Possible values are \code{"perc"}, \code{"norm"} and 
    \code{"basic"} (see \code{\link[boot]{boot.ci}}).}
  \item{alpha}{a numeric value giving the significance level to be used for 
    computing the confidence interval(s) (i.e., the confidence level is 
    \eqn{1 - }\code{alpha}), or \code{NULL}.}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
  \item{\dots}{if \code{bootType} is \code{"calibrate"}, additional arguments 
    to be passed to \code{\link{calibWeights}}.}
}
\value{
  An object of the same class as \code{indicator} is returned.  See 
  \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}} or 
  \code{\link{gini}} for details on the components.
}
\note{
  This function gives reasonable variance estimates for basic sample designs 
  such as simple random sampling or stratified simple random sampling.
}
\author{Andreas Alfons}
\seealso{
\code{\link{variance}}, \code{\link{calibWeights}}, \code{\link{arpr}}, 
\code{\link{qsr}}, \code{\link{rmpg}}, \code{\link{gini}}
}
\examples{
data(eusilc)
a <- arpr("eqIncome", weights = "rb050", data = eusilc)

## naive bootstrap
bootVar("eqIncome", weights = "rb050", design = "db040", 
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
bootVar("eqIncome", weights = "rb050", design = "db040", 
    data = eusilc, indicator = a, X = aux, totals = totals, 
    seed = 123)
}
\keyword{survey}
