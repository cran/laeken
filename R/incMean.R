# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# income mean
incMean <- function(inc, weights = NULL, years = NULL, 
        data = NULL, na.rm = FALSE) {
	## initializations
    if(!is.null(data)) {
        inc <- data[, inc]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(years)) years <- data[, years]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
    n <- length(inc)
    if(!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if(!is.null(years) && !is.numeric(years)) {
        stop("'years' must be a numeric vector")
    }
    if(is.null(data)) {  # check vector lengths
        if(!is.null(weights) && length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        if(!is.null(years) && length(years) != n) {
            stop("'years' must have the same length as 'x'")
        }
    }
    ## computations
    if(is.null(years)) {  # no breakdown
        xn <- weightedMean(inc, weights, na.rm=na.rm)
    } else {  # breakdown by years
        # define wrapper functions
        calcMean <- function(y, inc, weights, years, na.rm) {
            i <- years == y
            weightedMean(inc[i], weights[i], na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        xn <- sapply(ys, calcMean, inc=inc, 
            weights=weights, years=years, na.rm=na.rm)
        names(xn) <- ys  # use years as names
    }
    ## return results
    return(xn)
}
