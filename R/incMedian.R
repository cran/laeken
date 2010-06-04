# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# income median
incMedian <- function(inc, weights = NULL, sort = NULL, 
        years = NULL, data = NULL, na.rm = FALSE) {
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
    ## sort values and weights
    order <- if(is.null(sort)) order(inc) else order(inc, sort)
    inc <- inc[order]
    weights <- weights[order]  # also works if 'weights' is NULL
    ## computations
    if(is.null(years)) {  # no breakdown
        med <- weightedMedian(inc, weights, sorted=TRUE, na.rm=na.rm)
    } else {  # breakdown by years
        years <- years[order]
        # define wrapper functions
        calcMedian <- function(y, inc, weights, years, na.rm) {
            i <- years == y
            weightedMedian(inc[i], weights[i], sorted=TRUE, na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        med <- sapply(ys, calcMedian, inc=inc, 
            weights=weights, years=years, na.rm=na.rm)
        names(med) <- ys  # use years as names
    }
    ## return results
    return(med)
}
