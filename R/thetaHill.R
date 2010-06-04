# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

thetaHill <- function (x, k) {
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
    ## return Hill estimate
    k/(sum(log(x[(n-k+1):n])) - k*log(x0))
}
