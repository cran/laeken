\name{incMean}
\alias{incMean}
\title{Weighted mean income}
\description{
  Compute the weighted mean income.
}
\usage{
incMean(inc, weights = NULL, years = NULL, 
        data = NULL, na.rm = FALSE)
}
\arguments{
  \item{inc}{either a numeric vector giving the (equivalized disposable) income, 
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
  \item{data}{an optional \code{data.frame}.}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
}
\value{
  A numeric vector containing the value(s) of the weighted mean income is 
  returned. 
}
\author{Andreas Alfons}
\seealso{\code{\link{weightedMean}}}
\examples{
data(eusilc)
incMean("eqIncome", weights = "rb050", data = eusilc)
}
\keyword{survey}
