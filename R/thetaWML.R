# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#thetaWML <- function(x, k = NULL, x0 = NULL, 
#        weight = c("residuals", "probability"), const, 
#        bias = TRUE, tol = .Machine$double.eps^0.25, ...) {
#    
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
#    if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
#    x <- sort(x)  # sort values
#    n <- length(x)  # number of observations
#    if(haveK) {  # 'k' is supplied, threshold is determined
#        if(k >= n) stop("'k' must be smaller than the number of observed values")
#        x0 <- x[n-k]  # threshold (scale parameter)
#    } else {  # 'k' is not supplied, it is determined using threshold
#        # values are already sorted
#        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
#        k <- length(which(x > x0))
#    }
#    xt <- x[(n-k+1):n]  # tail (values larger than threshold)
#    y <- log(xt/x0)  # relative excesses
#    weight <- match.arg(weight)  # check type of robustness weights
#    
#    ## define robustness weight function and function for root finding
#    ## derivative of log(f) with respect to theta: 1/theta - log(xt/x0)
#    if(weight == "residuals") {
#        ## check tuning constant
#        if(missing(const)) const <- 2.5
#        else if(!is.numeric(const) || length(const) == 0) {
#            stop("'const' must be a numeric value")
#        } else const <- const[1]
#        ## some temporary values
#        h <- k:1
#        hy <- log(h/(k+1))
#        hsig <- sqrt(cumsum(1/h^2))
#        ## objective function
#        zeroTheta <- function(theta) {
#            r <- (theta*y + hy) / hsig  # standardized residuals
#            u <- pmin(1, const/abs(r))  # robustness weights
#            dlogf <- 1/theta - y  # derivative of log(f)
#            sum(u * dlogf)
#        }
#    } else {
#        ## check tuning constants
#        if(missing(const)) const <- rep.int(0.005, 2)
#        else if(!is.numeric(const) || length(const) == 0) {
#            stop("'const' must be a numeric vector of length two")
#        } else const <- rep(const, length.out=2)
#        p1 <- const[1]
#        p2 <- const[2]
#        ## objective function
#        zeroTheta <- function(theta) {
#            F <- 1 - (xt/x0)^(-theta)  # distribution function
#            u <- ifelse(F < p1, F/p1, ifelse(F <= 1-p2, 1, (1-F)/p2))  # robustness weights
#            dlogf <- 1/theta - y  # derivative of log(f)
#            sum(u * dlogf)
#        }
#    }
#    
#    ## solving sum(phi(xt,theta))=0
#    localUniroot <- function(f, interval = NULL, tol, ...) {
#        if(is.null(interval)) {
#            p <- if(haveK) thetaHill(x, k) else thetaHill(x, x0=x0)
#            interval <- c(0 + tol, 5 * p)  # default interval
#        }
#        uniroot(f, interval, ...)
#    }
#    theta <- localUniroot(zeroTheta, tol=tol, ...)$root
#    
#    ## optional bias correction
#    if(bias) {
#        if(weight == "residuals") {
#            r <- (theta*y + hy) / hsig  # standardized residuals
#            u <- pmin(1, const/abs(r))  # robustness weights
#            F <- 1 - (xt/x0)^(-theta)  # distribution function
#            deltaF <- diff(c(0, F))  # difference operator applied to F
#            dlogf <- 1/theta - y  # derivative of log(f)
#            d2logf <- -1/theta^2  # second derivative of log(f)
#            # derivative of robustness weight function
#            du <- ifelse(u == 1, 0, (-const)*y*hsig / (theta*y + hy)^2)
#            # bias correction term
#            bcorr <- -sum(u*dlogf*deltaF)/sum((du*dlogf + u*d2logf) * deltaF)
#        } else {
#            cp1 <- 1-p1
#            cp2 <- 1-p2
#            # bias correction term
#            bcorr <- (theta/2) * 
#                (2*cp1^2*log(cp1) + p1*cp1 + p1*cp2 + 2*p1*p2*log(p2)) / 
#                ((cp1*log(cp1))^2 - p1*cp1 - p1*cp2 + p1*p2*(log(p2))^2)
#        }
#        ## apply bias correction to theta
#        theta <- theta - bcorr
#    }
#    
#    ## return WML-estimate
#    theta
#}


thetaWML <- function(x, k = NULL, x0 = NULL, 
        weight = c("residuals", "probability"), 
        const, bias = TRUE, ...) {
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
    x <- sort(x)  # sort values
    if(missing(const)) .thetaWML(x, k, x0, weight, bias=bias, ...)
    else .thetaWML(x, k, x0, weight, const, bias, ...)
}

# internal function that assumes that data are ok and sorted
.thetaWML <- function(x, k = NULL, x0 = NULL, 
        weight = c("residuals", "probability"), const, 
        bias = TRUE, tol = .Machine$double.eps^0.25, ...) {
    
    n <- length(x)  # number of observations
    haveK <- !is.null(k)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- x[n-k]  # threshold (scale parameter)
    } else {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    xt <- x[(n-k+1):n]  # tail (values larger than threshold)
    y <- log(xt/x0)  # relative excesses
    weight <- match.arg(weight)  # check type of robustness weights
    
    ## define robustness weight function and function for root finding
    ## derivative of log(f) with respect to theta: 1/theta - log(xt/x0)
    if(weight == "residuals") {
        ## check tuning constant
        if(missing(const)) const <- 2.5
        else if(!is.numeric(const) || length(const) == 0) {
            stop("'const' must be a numeric value")
        } else const <- const[1]
        ## some temporary values
        h <- k:1
        hy <- log(h/(k+1))
        hsig <- sqrt(cumsum(1/h^2))
        ## objective function
        zeroTheta <- function(theta) {
            r <- (theta*y + hy) / hsig  # standardized residuals
            u <- pmin(1, const/abs(r))  # robustness weights
            dlogf <- 1/theta - y  # derivative of log(f)
            sum(u * dlogf)
        }
    } else {
        ## check tuning constants
        if(missing(const)) const <- rep.int(0.005, 2)
        else if(!is.numeric(const) || length(const) == 0) {
            stop("'const' must be a numeric vector of length two")
        } else const <- rep(const, length.out=2)
        p1 <- const[1]
        p2 <- const[2]
        ## objective function
        zeroTheta <- function(theta) {
            F <- 1 - (xt/x0)^(-theta)  # distribution function
            u <- ifelse(F < p1, F/p1, ifelse(F <= 1-p2, 1, (1-F)/p2))  # robustness weights
            dlogf <- 1/theta - y  # derivative of log(f)
            sum(u * dlogf)
        }
    }
    
    ## solving sum(phi(xt,theta))=0
    localUniroot <- function(f, interval = NULL, tol, ...) {
        if(is.null(interval)) {
            p <- if(haveK) .thetaHill(x, k) else .thetaHill(x, x0=x0)
            interval <- c(0 + tol, 5 * p)  # default interval
        }
        uniroot(f, interval, ...)
    }
    theta <- localUniroot(zeroTheta, tol=tol, ...)$root
    
    ## optional bias correction
    if(bias) {
        if(weight == "residuals") {
            r <- (theta*y + hy) / hsig  # standardized residuals
            u <- pmin(1, const/abs(r))  # robustness weights
            F <- 1 - (xt/x0)^(-theta)  # distribution function
            deltaF <- diff(c(0, F))  # difference operator applied to F
            dlogf <- 1/theta - y  # derivative of log(f)
            d2logf <- -1/theta^2  # second derivative of log(f)
            # derivative of robustness weight function
            du <- ifelse(u == 1, 0, (-const)*y*hsig / (theta*y + hy)^2)
            # bias correction term
            bcorr <- -sum(u*dlogf*deltaF)/sum((du*dlogf + u*d2logf) * deltaF)
        } else {
            cp1 <- 1-p1
            cp2 <- 1-p2
            # bias correction term
            bcorr <- (theta/2) * 
                (2*cp1^2*log(cp1) + p1*cp1 + p1*cp2 + 2*p1*p2*log(p2)) / 
                ((cp1*log(cp1))^2 - p1*cp1 - p1*cp2 + p1*p2*(log(p2))^2)
        }
        ## apply bias correction to theta
        theta <- theta - bcorr
    }
    
    ## return WML-estimate
    theta
}
