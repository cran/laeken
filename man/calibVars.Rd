\name{calibVars}
\alias{calibVars}
\docType{data}
\title{
Construct a matrix of binary variables for calibration
}
\description{
Construct a matrix of binary variables for calibration of sample weights 
according to known marginal population totals.
}
\usage{
calibVars(x)
}
\arguments{
  \item{x}{a vector that can be interpreted as factor, or a matrix or 
    \code{data.frame} consisting of such variables.}
}
\value{
A matrix of binary variables that indicate membership to the corresponding 
factor levels.
}
\author{Andreas Alfons}
\seealso{
\code{\link{calibWeights}}
}
\examples{
data(eusilc)
# default method
aux <- calibVars(eusilc$rb090)
head(aux)
# data.frame method
aux <- calibVars(eusilc[, c("db040", "rb090")])
head(aux)
}
\keyword{survey}
