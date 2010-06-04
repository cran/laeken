# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

thetaMoment <- function(x, k) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(!is.numeric(k) || length(k) == 0 || k[1] < 1) {
        stop("'k' must be a positive integer")
    } else k <- k[1]
    if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
    x <- sort(x)
    n <- length(x)
    if(k >= n) stop("'k' must be smaller than the number of observed values")
    x0 <- x[n-k]  # threshold (scale parameter)
    y <- log(x[(n-k+1):n]/x0)  # relative excesses
    ## moments
    M1 <- sum(y)/k  # first moment
    M2 <- sum(y^2)/k  # second moment
    ## moment estimator
    1/(M1 + 1 - 1/(2*(1-M1^2/M2)))
}
