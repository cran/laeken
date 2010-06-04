\name{thetaPDC}
\alias{thetaPDC}
\title{Partial density component (PDC) estimator}
\description{
  The partial density component (PDC) estimator estimates the shape parameter 
  of a Pareto distribution based on the relative excesses of observations 
  above a certain threshold.
}
\usage{
thetaPDC(x, k, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{the number of observations in the upper tail to which the Pareto 
    distribution is fitted.}
  \item{\dots}{additional arguments to be passed to \code{\link[stats]{nlm}} 
    (see \dQuote{Details}).}
}
\details{
  The PDC estimator is obtained by minimizing the integrated squared error 
  (ISE) criterion with an incomplete density mixture model.  The minimization 
  is carried out using \code{\link[stats]{nlm}}.  By default, the starting 
  value is the Hill estimator.
}
\value{
  The estimated shape parameter.
}
\references{
Vandewalle, B., Beirlant, J., Christmann, A., and Hubert, M. (2007) A robust 
estimator for the tail index of Pareto-type distributions. \emph{Computational 
Statistics & Data Analysis}, \bold{51}(12), 6252--6268.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}, \code{\link{thetaISE}}, \code{\link{thetaHill}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaPDC(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
