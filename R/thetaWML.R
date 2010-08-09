# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

thetaWML <- function(x, k, weight = c("residuals", "probability"), 
        const, bias = TRUE, tol = .Machine$double.eps^0.25, ...) {
    
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
    xt <- x[(n-k+1):n]  # tail (values larger than threshold)
    y <- log(xt/x0)  # relative excesses
    weight <- match.arg(weight)
    
    ## define weight function and function for root finding
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
            w <- pmin(1, const/abs(r))  # weights
            ## derivative of log(f) with respect to theta: 1/theta - log(xt/x0)
            sum(w * (1/theta - y))
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
            w <- ifelse(F < p1, F/p1, ifelse(F <= 1-p2, 1, (1-F)/p2))  # weights
            ## derivative of log(f) with respect to theta: 1/theta - log(xt/x0)
            sum(w * (1/theta - y))
        }
    }
    
#    ## solving sum(phi(xt,theta))=0
#    ## minimization of quadratic function instead of root finding
#    localNlm <- function(f, p = NULL, tol, ...) {
#        if(is.null(p)) p <- thetaHill(x, k)  # starting parameter
#        theta <- nlm(function(theta) f(theta)^2, p, ...)$estimate
#        if(abs(f(theta)) < tol) theta
#        else {
#            warning("M-estimator for 'theta' could not ", 
#                "be obtained, starting value is returned")
#            p
#        }
#    }
#    theta <- localNlm(zeroTheta, tol=tol, ...) 
    localUniroot <- function(f, interval = NULL, tol, ...) {
        if(is.null(interval)) {
            p <- thetaHill(x, k)  # Hill estimator
            interval <- c(0 + tol, 3 * p)  # default interval
        }
        uniroot(f, interval, ...)
    }
    theta <- localUniroot(zeroTheta, tol=tol, ...)$root
    
    ## optional bias correction
    if(bias) {
        if(weight == "residuals") {
            r <- (theta*y + hy) / hsig  # standardized residuals
            w <- pmin(1, const/abs(r))  # weights
            F <- 1 - (xt/x0)^(-theta)  # distribution function
            deltaF <- diff(c(0, F))  # difference operator applied to F
            dlogf <- 1/theta - y  # derivative of log(f)
            d2logf <- -1/theta^2  # second derivative of log(f)
            # derivative of weight function
            dw <- ifelse(w == 1, 0, (-const)*y*hsig / (theta*y + hy)^2)
            # bias correction
            bcorr <- -sum(w*dlogf*deltaF)/sum((dw*dlogf + w*d2logf) * deltaF)
    } else {
            cp1 <- 1-p1
            cp2 <- 1-p2
            # bias correction
            bcorr <- (theta/2) * 
                (2*cp1^2*log(cp1) + p1*cp1 + p1*cp2 + 2*p1*p2*log(p2)) / 
                ((cp1*log(cp1))^2 - p1*cp1 - p1*cp2 + p1*p2*(log(p2))^2)
        }
        ## apply bias correction to theta
        theta <- theta - bcorr
    }
    
    ## return WML estimate
    theta
}
    