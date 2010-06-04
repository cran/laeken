\name{fitPareto}
\alias{fitPareto}
\title{Fit income distribution models with the Pareto distribution}
\description{
  Fit a Pareto distribution to the upper tail of the income.  Since a 
  theoretical distribution is used for the upper tail, this is a 
  semi-parametric approach.
}
\usage{
fitPareto(x, k, method = "thetaPDC", groups = NULL, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{the number of observations in the upper tail to which the Pareto 
    distribution is fitted.}
  \item{method}{either a function or a character string specifying the function 
    to be used to estimate the shape parameter of the Pareto distibution, such 
    as \code{\link{thetaPDC}} (the default).  See \dQuote{Details} for 
    requirements for such a function and \dQuote{See also} for available 
    functions.}
  \item{groups}{a vector or factor specifying groups of elements of \code{x} 
    (e.g., households).  If supplied, each group of observations is expected to 
    have the same value in \code{x} (e.g., household income).  Only the values 
    of the first group member to appear are used for fitting the Pareto 
    distribution. For each group above the threshold, every group member is 
    assigned the same value.}
  \item{\dots}{addtional arguments to be passed to the specified method.}
}
\details{
  The function supplied to \code{method} should take a numeric vector (the 
  observations) as its first argument.  Its second argument should be the 
  number of observations in the upper tail to which the Pareto distribution is 
  fitted.  Additional arguments are passed via the \dots argument.
}
\value{
  A numeric vector with a Pareto distribution fit to the upper tail.
}
\author{Andreas Alfons and Josef Holzer}
\seealso{
  \code{\link{thetaPDC}}, \code{\link{thetaWML}}, \code{\link{thetaHill}}, 
  \code{\link{thetaISE}}, \code{\link{thetaLS}}, \code{\link{thetaMoment}}, 
  \code{\link{thetaQQ}}, \code{\link{thetaTM}}
}
\examples{
data(eusilc)

# gini coefficient without Pareto tail modeling
gini("eqIncome", weights = "rb050", data = eusilc)

# gini coefficient with Pareto tail modeling
eqIncome <- fitPareto(eusilc$eqIncome, 
    k = 75, groups = eusilc$db030)
gini(eqIncome, weights = eusilc$rb050)
}
\keyword{manip}
