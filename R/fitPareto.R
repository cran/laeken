# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

fitPareto <- function(x, k, method = "thetaPDC", groups = NULL, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(!is.numeric(k) || length(k) == 0 || k[1] < 1) {
        stop("'k' must be a positive integer")
    } else k <- k[1]
    if(is.null(groups)) values <- x
    else {
        if(!is.vector(groups) && !is.factor(groups)) {
            stop("'groups' must be a vector or factor")
        }
        if(length(groups) != length(x)) {
            stop("'groups' must have the same length as 'x'")
        }
        if(any(is.na(groups))) stop("'groups' contains missing values")
        unique <- !duplicated(groups)
        values <- x[unique]
    }
    ## check for missing values
    indices <- 1:length(values)
    if(any(i <- is.na(values))) indices <- indices[!i]
    ## order of observed values
    order <- order(values[indices])
    indicesSorted <- indices[order]  # indices of sorted vector
    n <- length(indicesSorted)
    if(k >= n) stop("'k' must be smaller than the number of observed values")
    xm <- values[indicesSorted[n-k]]  # threshold (scale parameter)
    ## estimate shape parameter
    theta <- do.call(method, list(values[indices], k, ...))
    ## fit Pareto distribution
    valuesPareto <- xm/runif(k)^(1/theta)
#    values[indicesSorted[(n-k+1):n]] <- valuesPareto
    values[indicesSorted[(n-k+1):n]] <- sort(valuesPareto)
    ## return values
    if(!is.null(groups)) {
        groups <- as.character(groups)
        names(values) <- groups[unique]
        values <- values[groups]
        names(values) <- names(x)
    }
    values
}
