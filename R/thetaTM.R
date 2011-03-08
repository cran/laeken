# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

thetaTM <- function(x, k = NULL, x0 = NULL, beta = 0.05) {
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
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- x[n-k]  # threshold (scale parameter)
    } else {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    # check trimming proportions
    if(length(beta) == 0) stop("'beta' must be a numeric vector of length two")
    else beta <- rep(beta, length.out=2)
    if(beta[1] < 0 || beta[1] >= 1) {
        stop("beta[1] (the trimming proportion for the lower end) ", 
            "must be greater or equal to 0 and smaller than 1")
    }
    if(beta[2] < 0 || beta[2] >= 1-beta[1]) {
        stop("beta[2] (the trimming proportion for the upper end) ", 
            "must be greater or equal to 0 and smaller than 1-beta[1] ",
            "(the trimming proportion for the lower end)")
    }
    # trimming
    kl <- trunc(k*beta[1])+1
    kh <- k - trunc(k*beta[2])
    i <- kl:kh
    c <- rep.int(0, k)
    c[i] <- 1/sum(cumsum(1/(k - 1:kh + 1))[i])
    # estimate theta
    1/sum(c*log(x[(n-k+1):n]/x0))
}
