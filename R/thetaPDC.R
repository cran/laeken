# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#thetaPDC <- function(x, k = NULL, x0 = NULL, w = NULL, 
#        tol = .Machine$double.eps^0.25, 
#        ...) {
#    ## initializations
#    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
#    haveK <- !is.null(k)
#    if(haveK) {  # if 'k' is supplied, it is always used
#        if(!is.numeric(k) || length(k) == 0 || k[1] < 1) {
#            stop("'k' must be a positive integer")
#        } else k <- k[1]
#    } else if(!is.null(x0)) {  # otherwise 'x0' (threshold) is used
#        if(!is.numeric(x0) || length(x0) == 0) stop("'x0' must be numeric")
#        else x0 <- x0[1]
#    } else stop("either 'k' or 'x0' must be supplied")
#    haveW <- !is.null(w)
#    if(haveW) {  # sample weights are supplied
#        if(!is.numeric(w) || length(w) != length(x)) {
#            stop("'w' must be numeric vector of the same length as 'x'")
#        }
#        if(any(w < 0)) stop("negative weights in 'w'")
#        if(any(i <- is.na(x))) {  # remove missing values
#            x <- x[!i]
#            w <- w[!i]
#        }
#        # sort values and sample weights
#        order <- order(x)
#        x <- x[order]
#        w <- w[order]
#    } else {  # no sample weights
#        if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
#        x <- sort(x)  # sort values
#    }
#    n <- length(x)  # number of observations
#    if(haveK) {  # 'k' is supplied, threshold is determined
#        if(k >= n) stop("'k' must be smaller than the number of observed values")
#        x0 <- x[n-k]  # threshold (scale parameter)
#    } else {  # 'k' is not supplied, it is determined using threshold
#        # values are already sorted
#        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
#        k <- length(which(x > x0))
#    }
#    ## computations
#    y <- x[(n-k+1):n]/x0  # relative excesses
#    if(haveW) {
#        wTail <- w[(n-k+1):n]
#        ## weighted integrated squared error distance criterion with incomplete 
#        ## density mixture model
#        # w ... sample weights
#        # u ... robustness weights (from incomplete density mixture model)
#        ISE <- function(theta, y, w) {
#            f <- theta*y^(-1-theta)
#            wm <- weighted.mean(f, w)  # weighted mean as unbiased estimator of expectation of f
#            pf2 <- theta^2/(2*theta+1)  # primitive of f^2
#            u <- wm/pf2
#            u^2*pf2 - 2*u*wm
#        }
#    } else {
#        wTail <- NULL
#        ## integrated squared error distance criterion with incomplete density 
#        ## mixture model
#        # w ... sample weights (not needed here, only available to have the 
#        #       same function definition)
#        # u ... robustness weights (from incomplete density mixture model)
#        ISE <- function(theta, y, w) {
#            f <- theta*y^(-1-theta)
#            m <- mean(f)  # mean as unbiased estimator of expectation of f
#            pf2 <- theta^2/(2*theta+1)  # primitive of f^2
#            u <- m/pf2
#            u^2*pf2 - 2*u*m
#        }
#    }
##    ## minimize
##    localNlm <- function(f, p = NULL, ...) {
##        if(is.null(p)) {  # obtain default starting parameter
##            p <- if(haveK) thetaHill(x, k, w=w) else thetaHill(x, x0=x0, w=w)
##        }
##        nlm(f, p, ...)
##    }
##    localNlm(ISE, y=y, w=wTail, ...)$estimate
#    ## optimize
#    localOptimize <- function(f, interval = NULL, tol, ...) {
#        if(is.null(interval)) {
#            p <- if(haveK) thetaHill(x, k, w=w) else thetaHill(x, x0=x0, w=w)
#            interval <- c(0 + tol, 3 * p)  # default interval
#        }
#        optimize(f, interval, ...)
#    }
#    localOptimize(ISE, y=y, w=wTail, tol=tol, ...)$minimum
#}


thetaPDC <- function(x, k = NULL, x0 = NULL, w = NULL, ...) {
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
    haveW <- !is.null(w)
    if(haveW) {  # sample weights are supplied
        if(!is.numeric(w) || length(w) != length(x)) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
        if(any(i <- is.na(x))) {  # remove missing values
            x <- x[!i]
            w <- w[!i]
        }
        # sort values and sample weights
        order <- order(x)
        x <- x[order]
        w <- w[order]
    } else {  # no sample weights
        if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
        x <- sort(x)  # sort values
    }
    .thetaPDC(x, k, x0, w, ...)
}

# internal function that assumes that data are ok and sorted
.thetaPDC <- function(x, k = NULL, x0 = NULL, w = NULL, 
        tol = .Machine$double.eps^0.25, ...) {
    n <- length(x)  # number of observations
    haveK <- !is.null(k)
    haveW <- !is.null(w)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- x[n-k]  # threshold (scale parameter)
    } else {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    ## computations
    y <- x[(n-k+1):n]/x0  # relative excesses
    if(haveW) {
        wTail <- w[(n-k+1):n]
        ## weighted integrated squared error distance criterion with incomplete 
        ## density mixture model
        # w ... sample weights
        # u ... robustness weights (from incomplete density mixture model)
        ISE <- function(theta, y, w) {
            f <- theta*y^(-1-theta)
            wm <- weighted.mean(f, w)  # weighted mean as unbiased estimator of expectation of f
            pf2 <- theta^2/(2*theta+1)  # primitive of f^2
            u <- wm/pf2
            u^2*pf2 - 2*u*wm
        }
    } else {
        wTail <- NULL
        ## integrated squared error distance criterion with incomplete density 
        ## mixture model
        # w ... sample weights (not needed here, only available to have the 
        #       same function definition)
        # u ... robustness weights (from incomplete density mixture model)
        ISE <- function(theta, y, w) {
            f <- theta*y^(-1-theta)
            m <- mean(f)  # mean as unbiased estimator of expectation of f
            pf2 <- theta^2/(2*theta+1)  # primitive of f^2
            u <- m/pf2
            u^2*pf2 - 2*u*m
        }
    }
    ## optimize
    localOptimize <- function(f, interval = NULL, tol, ...) {
        if(is.null(interval)) {
            p <- if(haveK) .thetaHill(x, k, w=w) else .thetaHill(x, x0=x0, w=w)
            interval <- c(0 + tol, 3 * p)  # default interval
        }
        optimize(f, interval, ...)
    }
    localOptimize(ISE, y=y, w=wTail, tol=tol, ...)$minimum
}
