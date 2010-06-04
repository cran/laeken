# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## gini coefficient
gini <- function(inc, weights = NULL, sort = NULL, years = NULL, 
        breakdown = NULL, design = NULL, data = NULL, var = NULL, 
        alpha = 0.05, na.rm = FALSE, ...) {
    ## initializations
    byYear <- !is.null(years)
    byStratum <- !is.null(breakdown)
    if(!is.null(data)) {
        inc <- data[, inc]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(sort)) sort <- data[, sort]
        if(byYear) years <- data[, years]
        if(byStratum) breakdown <- data[, breakdown]
        if(!is.null(var) && !is.null(design)) design <- data[, design]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
    n <- length(inc)
    if(is.null(weights)) weights <- weights <- rep.int(1, n)
    else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    if(!is.null(sort) && !is.vector(sort) && !is.ordered(sort)) {
        stop("'sort' must be a vector or ordered factor")
    }
    if(byYear && !is.numeric(years)) {
        stop("'years' must be a numeric vector")
    }
    if(byStratum) {
        if(!is.vector(breakdown) && !is.factor(breakdown)) {
            stop("'breakdown' must be a vector or factor")
        } else breakdown <- as.factor(breakdown)
    }
    if(is.null(data)) {  # check vector lengths
        if(length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        if(!is.null(sort) && length(sort) != n) {
            stop("'sort' must have the same length as 'x'")
        }
        if(byYear && length(years) != n) {
            stop("'years' must have the same length as 'x'")
        }
        if(byStratum && length(breakdown) != n) {
            stop("'breakdown' must have the same length as 'x'")
        }
    }
    ## computations
    # Gini by year (if requested)
    if(byYear) {
        ys <- sort(unique(years))  # unique years
        gc <- function(y, inc, weights, sort, years, na.rm) {
            i <- years == y
            giniCoeff(inc[i], weights[i], sort[i], na.rm=na.rm)
        }
        value <- sapply(ys, gc, inc=inc, weights=weights, 
            sort=sort, years=years, na.rm=na.rm)
        names(value) <- ys  # use years as names
    } else {
        ys <- NULL
        value <- giniCoeff(inc, weights, sort, na.rm=na.rm)
    }
    # Gini by stratum (if requested)
    if(byStratum) {
        gcR <- function(i, inc, weights, sort, na.rm) {
            giniCoeff(inc[i], weights[i], sort[i], na.rm=na.rm)
        }
        valueByStratum <- aggregate(1:n, 
            if(byYear) list(year=years, stratum=breakdown) else list(stratum=breakdown), 
            gcR, inc=inc, weights=weights, sort=sort, na.rm=na.rm)
        names(valueByStratum)[ncol(valueByStratum)] <- "value"
        rs <- levels(breakdown)  # unique strata
    } else valueByStratum <- rs <- NULL
    ## create object of class "qsr"
    res <- constructGini(value=value, 
        valueByStratum=valueByStratum, 
        years=ys, strata=rs)
    # variance estimation (if requested)
    if(!is.null(var)) {
        res <- variance(inc, weights, years, breakdown, design, 
            indicator=res, alpha=alpha, na.rm=na.rm, type=var, ...)
    }
    ## return result
    return(res)
}

## workhorse
giniCoeff <- function(x, weights = NULL, sort = NULL, na.rm = FALSE) {
    # initializations
    if(isTRUE(na.rm)){
        indices <- !is.na(x)
        x <- x[indices]
        if(!is.null(weights)) weights <- weights[indices]
        if(!is.null(sort)) sort <- sort[indices]
    } else if(any(is.na(x))) return(NA)
    # sort values and weights
    order <- if(is.null(sort)) order(x) else order(x, sort)
    x <- x[order]  # order values
    if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
    else weights <- weights[order]  # order weights
    ## calculations
    wx <- weights * x  # weighted values
    sw <- sum(weights)  # sum of weights
    cw <- cumsum(weights)  # cumulative sum of weights
    100 * ((2 * sum(wx*cw) - sum(weights^2 * x)) / (sw * sum(wx)) - 1)
}
