# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

thetaQQ <- function(x, k = NULL, x0 = NULL) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    haveK <- !is.null(k)
    if(haveK) {  # if 'k' is supplied, it is always used
        if(!is.numeric(k) || length(k) == 0 || k[1] < 1) {
            stop("'k' must be a positive integer")
        } else k <- k[1]
    } else if(!is.null(x0)) {  # otherwise 'x0' (threshold) is used
        if(!is.numeric(x0) || length(x0) == 0) stop("'x0' must be numeric")
        else x0 <- x0[1]
    } else stop("either 'k' or 'x0' must be supplied")
    if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
    x <- sort(x)
    n <- length(x)
#    if(haveK) {  # 'k' is supplied, threshold is determined
#        if(k >= n) stop("'k' must be smaller than the number of observed values")
#        x0 <- x[n-k]  # threshold (scale parameter)
#    } else {  # 'k' is not supplied, it is determined using threshold
#        # values are already sorted
#        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
#        k <- length(which(x > x0))
#    }
    if(!haveK) {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    ## calculations
    logx <- log(x[n:(n-k+1)])  # lograrithm of reversed tail
    h <- -log((1:k)/(k+1))
    ## QQ-estimator
    (k*sum(h^2) - sum(h)^2) / sum(h * (k*logx - sum(logx)))
}
