# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# income quintile
incQuintile <- function(inc, weights = NULL, sort = NULL, 
        years = NULL, k = c(1, 4), data = NULL, na.rm = FALSE) {
    ## initializations
    if(!is.null(data)) {
        inc <- data[, inc]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(sort)) sort <- data[, sort]
        if(!is.null(years)) years <- data[, years]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
    n <- length(inc)
    if(!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if(!is.null(sort) && !is.vector(sort) && !is.ordered(sort)) {
        stop("'sort' must be a vector or ordered factor")
    }
    if(!is.null(years) && !is.numeric(years)) {
        stop("'years' must be a numeric vector")
    }
    if(is.null(data)) {  # check vector lengths
        if(!is.null(weights) && length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        if(!is.null(sort) && length(sort) != n) {
            stop("'sort' must have the same length as 'x'")
        }
        if(!is.null(years) && length(years) != n) {
            stop("'years' must have the same length as 'x'")
        }
    }
    if(!is.numeric(k) || length(k) == 0 || any(k < -0.5 | k >= 5.5)) {
        stop("'k' must be a vector of integers between 0 and 5")
    } else k <- round(k)
    ## sort values and weights
    order <- if(is.null(sort)) order(inc) else order(inc, sort)
    inc <- inc[order]
    weights <- weights[order]  # also works if 'weights' is NULL
    ## computations
    if(is.null(years)) {  # no breakdown
        q <- weightedQuantile(inc, weights, probs=k/5, sorted=TRUE, na.rm=na.rm)
        names(q) <- k  # use quintile numbers as names
    } else {  # breakdown by years
        years <- years[order]
        # define wrapper functions
        calcQuantile <- function(y, inc, weights, years, k, na.rm) {
            i <- years == y
            weightedQuantile(inc[i], weights[i], 
                probs=k/5, sorted=TRUE, na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        q <- t(sapply(ys, calcQuantile, inc=inc, 
                weights=weights, years=years, k=k, na.rm=na.rm))
        rownames(q) <- ys  # use years as row names
        colnames(q) <- k   # use quintile numbers as column names
    }
    ## return results
    return(q)
}
