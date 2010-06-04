# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

thetaQQ <- function(x, k) {
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
    logx <- log(x[n:(n-k+1)])  # lograrithm of reversed tail
    h <- -log((1:k)/(k+1))
    ## QQ-estimator
    (k*sum(h^2) - sum(h)^2) / sum(h * (k*logx - sum(logx)))
}
