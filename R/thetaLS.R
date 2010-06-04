# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

## should we return estimate for x0? if so, don't we need to re-estimate theta?
## => iterative procedure until change smaller than a threshold?

thetaLS <- function(x, k) {
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
    z <- log(x[(n-k+1):n])
    zm <- mean(z)
    pk <- c((1:(k-1))/k, k/(k+1))  # regression parameters
    ck <- -log(1-pk)
    ckm <- mean(ck)
    ## LS estimator
    mean((ck - ckm)^2) / (mean(ck*z) - ckm*zm)
}
