\name{thetaWML}
\alias{thetaWML}
\title{Weighted maximum likelihood estimator}
\description{
  Estimate the shape parameter of a Pareto distribution using a weighted 
  maximum likelihood approach.
}
\usage{
thetaWML(x, k, weight = c("residuals", "probability"), const, 
         bias = TRUE, tol = .Machine$double.eps^0.25, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{the number of observations in the upper tail to which the Pareto 
    distribution is fitted.}
  \item{weight}{a character string specifying the weight function to be used.  
    If \code{"residuals"} (the default), the weight function is based on 
    standardized residuals.  If \code{"probability"}, probability based 
    weighting is used.  Partial string matching allows these names to be 
    abbreviated.}
  \item{const}{Tuning constant(s) that control the robustness of the method.  
    If \code{weight="residuals"}, a single numeric value is required (the 
    default is 2.5).  If \code{weight="probability"}, a numeric vector of 
    length two must be supplied (a single numeric value is recycled; the 
    default is 0.005 for both tuning parameters).  See the reference for more 
    details.}
  \item{bias}{a logical indicating whether bias correction should be applied.}
  \item{tol}{the desired accuracy.}
  \item{\dots}{additional arguments to be passed to 
    \code{\link[stats]{uniroot}} (see \dQuote{Details}).}
}
\details{
  The weighted maximum likelihood estimator belongs to the class of 
  M-estimators.  In order to obtain the estimator, the root of a certain 
  function needs to be found, which is implemented using 
  \code{\link[stats]{uniroot}}.
}
\value{
  The estimated shape parameter.
}
\references{
Dupuis, D.J. and Morgenthaler, S. (2002) Robust weighted likelihood estimators 
with an application to bivariate extreme value problems. \emph{The Canadian 
Journal of Statistics}, \bold{30}(1), 17--36.

Dupuis, D.J. and Victoria-Feser, M.-P. (2006) A robust prediction error 
criterion for Pareto modelling of upper tails. \emph{The Canadian Journal of 
Statistics}, \bold{34}(4), 639--658.
}
\author{Josef Holzer and Andreas Alfons}
\seealso{
  \code{\link{fitPareto}}
}
\examples{
data(eusilc)
# equivalized disposable income is equal for each household
# member, therefore only one household member is taken
thetaWML(eusilc$eqIncome[!duplicated(eusilc$db030)], k = 75)
}
\keyword{manip}
