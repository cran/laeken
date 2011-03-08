# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

fitPareto <- function(x, k = NULL, x0 = NULL, 
        method = "thetaPDC", groups = NULL, w = NULL, ...) {
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
    nam <- argNames(method)
    useW <- !is.null(w) && ("w" %in% nam)
    if(useW && (!is.numeric(w) || length(w) != length(x))) {
        stop("'w' must be numeric vector of the same length as 'x'")
    }
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
        values <- x[unique]
        if(useW) w <- w[unique]
    } else values <- x
    ## check for missing values
    indices <- 1:length(values)
    if(any(i <- is.na(values))) indices <- indices[!i]
    ## order of observed values
    order <- order(values[indices])
    indicesSorted <- indices[order]  # indices of sorted vector
    n <- length(indicesSorted)
    ## start constructing call to 'method' for estimation of shape parameter
    dots <- list(values[indices], ...)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- values[indicesSorted[n-k]]  # threshold (scale parameter)
        dots$k <- k  # 'method' is expected to have 'k' as argument
    } else {  # 'k' is not supplied, it is determined using threshold
        if(x0 >= values[indicesSorted[n]]) {  # compare to sorted values
            stop("'x0' must be smaller than the largest value")
        }
        k <- length(which(values[indices] > x0))  # number of observations in tail
        dots$x0 <- x0  # 'method' is expected to have threshold 'x0' as argument
    }
    ## estimate shape parameter
    if(useW) dots$w <- w[indices]
    theta <- do.call(method, dots)
    ## fit Pareto distribution
    valuesPareto <- x0/runif(k)^(1/theta)
    values[indicesSorted[(n-k+1):n]] <- sort(valuesPareto)
    ## return values
    if(haveGroups) {
        groups <- as.character(groups)
        names(values) <- groups[unique]
        values <- values[groups]
        names(values) <- names(x)
    }
    values
}
