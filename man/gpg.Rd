\name{gpg}
\alias{gpg}
\title{
Gender pay (wage) gap.
}
\description{
Estimate the gender pay (wage) gap.
}
\usage{
gpg(inc, gender = NULL, method = c("mean", "median"), 
    weights = NULL, sort = NULL, years = NULL, breakdown = NULL, 
    design = NULL, data = NULL, var = NULL, alpha = 0.05, 
    na.rm = FALSE, ...)
}
\arguments{
  \item{inc}{either a numeric vector giving the equivalized disposable income, 
    or (if \code{data} is not \code{NULL}) a character string, an integer or a 
    logical vector specifying the corresponding column of \code{data}.}
  \item{gender}{either a factor giving the gender, 
    or (if \code{data} is not \code{NULL}) a character string, an integer or a 
    logical vector specifying the corresponding column of \code{data}.}
  \item{method}{a character string specifying the method to be used.  Possible 
    values are \code{"mean"} for the mean, and \code{"median"} for the median.
    If weights are provided, the weighted mean or weighted median is estimated.}
  \item{weights}{optional; either a numeric vector giving the personal sample 
    weights, or (if \code{data} is not \code{NULL}) a character string, an 
    integer or a logical vector specifying the corresponding column of 
    \code{data}.}
  \item{sort}{optional; either a numeric vector giving the personal IDs to be 
    used as tie-breakers for sorting, or (if \code{data} is not \code{NULL}) a 
    character string, an integer or a logical vector specifying the 
    corresponding column of \code{data}.}
  \item{years}{optional; either a numeric vector giving the different years of 
    the survey, or (if \code{data} is not \code{NULL}) a character string, 
    an integer or a logical vector specifying the corresponding column of 
    \code{data}.  If supplied, values are computed for each year.}
  \item{breakdown}{optional; either a numeric vector giving different strata, 
    or (if \code{data} is not \code{NULL}) a character string, an integer or a 
    logical vector specifying the corresponding column of \code{data}.  If 
    supplied, the values for each stratum are computed in addition to the 
    overall value.}
  \item{design}{optional and only used if \code{var} is not \code{NULL}; either 
    an integer vector or factor giving different strata for stratified sampling 
    designs, or (if \code{data} is not \code{NULL}) a character string, an 
    integer or a logical vector specifying the corresponding column of 
    \code{data}.}
  \item{data}{an optional \code{data.frame}.}
  \item{var}{a character string specifying the type of variance estimation to 
    be used, or \code{NULL} to omit variance estimation.  See 
    \code{\link{variance}} for possible values.}
  \item{alpha}{numeric; if \code{var} is not \code{NULL}, this gives the 
    significance level to be used for computing the confidence interval 
    (i.e., the confidence level is \eqn{1 - }\code{alpha}).}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
  \item{\dots}{if \code{var} is not \code{NULL}, additional arguments to be 
    passed to \code{\link{variance}}.}
}
\details{
  The implementation strictly follows the Eurostat definition (with default 
  method \code{"mean"} and alternative method \code{"median"}).  If weights are 
  provided, the weighted mean or weighted median is estimated. 
}
\value{
  A list of class \code{"gpg"} (which inherits from the class 
  \code{"indicator"}) with the following components:
  \item{value}{a numeric vector containing the overall value(s).}
  \item{valueByStratum}{a \code{data.frame} containing the values by stratum, 
    or \code{NULL}.}
  \item{varMethod}{a character string specifying the type of variance 
    estimation used, or \code{NULL} if variance estimation was omitted.}
  \item{var}{a numeric vector containing the variance estimate(s), or 
    \code{NULL}.}
  \item{varByStratum}{a \code{data.frame} containing the variance estimates by 
    stratum, or \code{NULL}.}
  \item{ci}{a numeric vector or matrix containing the lower and upper endpoints 
    of the confidence interval(s), or \code{NULL}.}
   \item{ciByStratum}{a \code{data.frame} containing the lower and upper 
    endpoints of the confidence intervals by stratum, or \code{NULL}.}
  \item{alpha}{a numeric value giving the significance level used for computing 
    the confidence interv
    al(s) (i.e., the confidence level is 
    \eqn{1 - }\code{alpha}), or \code{NULL}.}
  \item{years}{a numeric vector containing the different years of the survey.}
  \item{strata}{a character vector containing the different strata of the 
    breakdown.}
}
\references{
Working group on Statistics on Income and Living Conditions (2004) Common 
cross-sectional EU indicators based on EU-SILC; the gender pay gap. 
\emph{EU-SILC 131-rev/04}, Eurostat.
}
\author{
Matthias Templ and Alexander Haider, using code for breaking down estimation by 
Andreas Alfons.
}
\seealso{
\code{\link{variance}}, \code{\link{qsr}}, \code{\link{gini}}
}
\examples{
data(eusilc)
### clearly, children and elder people may not work in Austria:
eusilc <- eusilc[eusilc$age > 17 & eusilc$age < 66, ]
### full time workers:
eusilc <- eusilc[eusilc$pl030 == 1, ]
### employees's cash income:
### py010n
eusilc <- eusilc[!is.na(eusilc$py010n), ]
### for estimation of the GPG, use hourly rates of people 
### who earn money and NOT just yearly income as done in 
### the following examples!

median_no_breakdown <- gpg("py010n", "rb090", method = "median", 
    weights = "rb050", data = eusilc)
mean_no_breakdown <- gpg("py010n", "rb090", method = "mean", 
    weights = "rb050", data = eusilc)

variance("py010n", gender = "rb090", method = "median", 
    weights = "rb050", design = "db040", data = eusilc, 
    indicator = median_no_breakdown, bootType = "naive", 
    seed = 123)

variance("py010n", gender = "rb090", method = "mean", 
    weights = "rb050", design = "db040", data = eusilc, 
    indicator = mean_no_breakdown, bootType = "naive", 
    seed = 123)


median_breakdown_area <- gpg("py010n", "rb090", 
    method = "median", breakdown = "db040", 
    weights = "rb050", data = eusilc)

variance("py010n", gender = "rb090", method = "median", 
    breakdown = "db040", weights = "rb050", design = "db040", 
    data = eusilc, indicator = median_breakdown_area, 
    bootType = "naive", seed = 123)


mean_breakdown_area <- gpg("py010n", "rb090", method = "mean", 
    breakdown = "db040", weights = "rb050", data = eusilc)

variance("py010n", gender = "rb090", method = "mean", 
    breakdown = "db040", weights = "rb050", design = "db040", 
    data = eusilc, indicator = mean_breakdown_area, 
    bootType = "naive", seed = 123)
}
\keyword{survey}
