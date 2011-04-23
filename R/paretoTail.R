# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

paretoTail <- function(x, k = NULL, x0 = NULL, method = "thetaPDC", 
        groups = NULL, w = NULL, alpha = 0.01, ...) {
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
    if(is.character(method)) method <- getDotTheta(method)
    nam <- argNames(method)
    useW <- !is.null(w) && ("w" %in% nam)
    if(useW && (!is.numeric(w) || length(w) != length(x))) {
        stop("'w' must be numeric vector of the same length as 'x'")
    }
    ## allow for cluster effect
    haveGroups <- !is.null(groups)
    if(haveGroups) {
        if(!is.vector(groups) && !is.factor(groups)) {
            stop("'groups' must be a vector or factor")
        }
        if(length(groups) != length(x)) {
            stop("'groups' must have the same length as 'x'")
        }
        if(any(is.na(groups))) stop("'groups' contains missing values")
        unique <- !duplicated(groups)
        xx <- x[unique]
        if(useW) ww <- w[unique]
    } else {
        xx <- x
        if(useW) ww <- w
    }
    xx <- unname(xx)
    ## check for missing values
    if(any(i <- is.na(xx))) {
        xx <- xx[!i]
        if(useW) ww <- ww[!i]
    }
    ## order of observed values
    order <- order(xx)
    xx <- xx[order]
    if(useW) ww <- ww[order]
    n <- length(xx)
    ## start constructing call to 'method' for estimation of shape parameter
    dots <- list(xx, ...)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- xx[n-k]  # threshold (scale parameter)
        dots$k <- k  # 'method' is expected to have 'k' as argument
    } else {  # 'k' is not supplied, it is determined using threshold
        if(x0 >= xx[n]) {  # compare to sorted values
            stop("'x0' must be smaller than the largest value")
        }
        k <- length(which(xx > x0))  # number of observations in tail
        dots$x0 <- x0  # 'method' is expected to have threshold 'x0' as argument
    }
    ## estimate shape parameter
    if(useW) dots$w <- ww
    theta <- do.call(method, dots)
    ## indicate observations in tail
    if(haveGroups) {
        tail <- groups[unique]
        tail <- tail[!i]
        tail <- tail[order]
        tail <- tail[(n-k+1):n]
    } else tail <- order[(n-k+1):n]
    ## flag suspicious observations (nonrepresentative outliers)
    if(!is.numeric(alpha) || length(alpha) == 0 || alpha < 0 || alpha > 1) {
        stop("'alpha' must be a numeric value in [0,1]")
    } else alpha <- alpha[1]
    q <- qpareto(1-alpha, x0, theta)  # quantile of the Pareto distribution
    if(haveGroups) {
        out <- which(xx[(n-k+1):n] > q)
        out <- tail[out]
    } else {
        out <- unname(which(x > q))
        out <- out[order(x[out])]
    }
    
    ## return object
    res <- list(x=x, k=k, groups=groups, w=w, method=method, 
        x0=x0, theta=theta, tail=tail, alpha=alpha, out=out)
    class(res) <- "paretoTail"
    res
}

replaceTail <- function(x, ...) UseMethod("replaceTail")

replaceTail.paretoTail <- function(x, all = TRUE, ...) {
    which <- if(isTRUE(all)) x$tail else x$out
    k <- length(which)  # number of observations to be replaced
    res <- x$x
    if(k > 0) {
        new <- sort(rpareto(k, x$x0, x$theta))
        groups <- x$groups
        if(is.null(groups)) res[which] <- new
        else {
            groups <- as.character(groups)
            which <- as.character(which)
            replace <- which(groups %in% which)
            names(new) <- which
            new <- new[groups[replace]]
            names(new) <- names(res[replace])
            res[replace] <- new
        }
    }
    res
}

#replaceOut <- function(x) UseMethod("replaceOut")
#
#replaceOut.paretoTail <- function(x) {
#    out <- x$out
#    nout <- length(out)  # number of nonrepresentative outliers
#    new <- sort(rpareto(nout, x$x0, x$theta))
#    res <- x$x
#    if(nout > 0) {
#        groups <- x$groups
#        if(is.null(groups)) res[out] <- new
#        else {
#            groups <- as.character(groups)
#            out <- as.character(out)
#            replace <- which(groups %in% out)
#            names(new) <- out
#            new <- new[groups[replace]]
#            names(new) <- names(res[replace])
#            res[replace] <- new
#        }
#    }
#    res
#}

#replaceOut <- function(x) replaceTail(x, all=FALSE)

replaceOut <- function(x, ...) {
    localReplaceTail <- function(x, all, ...) replaceTail(x, all=FALSE, ...)
    localReplaceTail(x, ...)
}


reweightOut <- function(x, ...) UseMethod("reweightOut")

reweightOut.paretoTail <- function(x, X, w = NULL, ...) {
    # in case of sample weights, set weights of outliers to one and calibrate 
    # other observations
    # otherwise, set weights of outliers to zero and weights of other 
    # observations to one
    out <- x$out
    n <- length(x$x)  # number of observations
    if(is.null(x$w)) {
        # check supplied weights
        if(!is.null(w) && (!is.numeric(w) || length(w) != n)) {
            stop(sprintf("'w' must be numeric vector of length %d", n))
        }
    } else w <- x$w
    if(length(out) > 0) {  # nonrepresentative outliers
        groups <- x$groups
        if(!is.null(groups)) out <- which(groups %in% out)
        if(is.null(w)) {
            w <- rep.int(1, n)
            w[out] <- 0
        } else {
            totals <- apply(X, 2, function(i) sum(i*w))
            args <- list(...)
            args$X <- X[-out, , drop=FALSE]
            args$d <- w[-out]
            w[out] <- 1  # set weight of nonrepresentative outliers to 1
            totalsOut <- apply(X[out, , drop=FALSE], 2, sum)
            args$totals <- totals - totalsOut
            g <- do.call("calibWeights", args)
            w[-out] <- g * args$d
        }
    }
    w
}


shrinkOut <- function(x, ...) UseMethod("shrinkOut")

shrinkOut.paretoTail <- function(x, ...) {
    # winsorize outliers in the upper tail
    out <- x$out
    res <- x$x
    if(length(out) > 0) {  # nonrepresentative outliers
        new <- qpareto(1-x$alpha, x$x0, x$theta)  # quantile of the Pareto distribution
        groups <- x$groups
        if(!is.null(groups)) out <- which(groups %in% out)
        res[out] <- new
    }
    res
}


## print method for class "paretoTail"
print.paretoTail <- function(x, ...) {
    cat("Threshold: ")
    cat(x$x0, ...)
    items <- if(is.null(x$groups)) "observations" else "groups"
    cat(sprintf("\nNumber of %s in the tail: ", items))
    cat(x$k, ...)
    cat("\nShape parameter: ")
    cat(x$theta, ...)
    cat(sprintf("\n\nOutlying %s:\n", items))
    print(x$out, ...)
}


## utility functions for Pareto distribution
dpareto <- function(x, x0 = 1, theta = 1) theta*x0^theta / x^(theta+1)
ppareto <- function(q, x0 = 1, theta = 1) 1 - (q/x0)^(-theta)
qpareto <- function(p, x0 = 1, theta = 1) unname(x0*exp(qexp(p)/theta))
rpareto <- function(n, x0 = 1, theta = 1) x0/runif(n)^(1/theta)


## other utility functions
getDotTheta <- function(method) {
    if(length(method) == 0) stop("'method' has length 0")
    else method <- method[1]
    if(method %in% c("thetaPDC", "thetaISE", "thetaWML", "thetaHill")) {
        method <- paste(".", method, sep="")
    }
    method
}
