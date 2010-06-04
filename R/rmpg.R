# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## relative median at-risk-of-poverty gap
rmpg <- function(inc, weights = NULL, sort = NULL, years = NULL, 
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
    if(byYear) {  # RMPG by year
        ys <- sort(unique(years))
        ts <- arpt(inc, weights, sort, years, na.rm=na.rm)  # thresholds
        rg <- function(y, t, inc, weights, sort, years, na.rm) {
            i <- years == y
            relativeGap(inc[i], weights[i], sort[i], t, na.rm=na.rm)
        }
        value <- mapply(rg, y=ys, t=ts, MoreArgs=list(inc=inc, 
                weights=weights, sort=sort, years=years, na.rm=na.rm))
        names(value) <- ys  # use years as names
        if(byStratum) {
            rg1 <- function(i, inc, weights, sort, years, ts, na.rm) {
                y <- years[i[1]]
                t <- ts[as.character(y)]
                relativeGap(inc[i], weights[i], sort[i], t, na.rm=na.rm)
            }
            valueByStratum <- aggregate(1:n, list(year=years, stratum=breakdown), 
                rg1, inc=inc, weights=weights, sort=sort, years=years, ts=ts, 
                na.rm=na.rm)
            names(valueByStratum)[3] <- "value"
        } else valueByStratum <- NULL
    } else {  # RMPG for only one year
        ys <- NULL
        ts <- arpt(inc, weights, sort, na.rm=na.rm)  # threshold
        value <- relativeGap(inc, weights, sort, ts, na.rm=na.rm)
        if(byStratum) {
            rg2 <- function(i, inc, weights, sort, ts, na.rm) {
                relativeGap(inc[i], weights[i], sort[i], ts, na.rm=na.rm)
            }
            valueByStratum <- aggregate(1:n, list(stratum=breakdown), 
                rg2, inc=inc, weights=weights, sort=sort, ts=ts, na.rm=na.rm)
            names(valueByStratum)[2] <- "value"
        } else valueByStratum <- NULL
    }
    rs <- levels(breakdown)  # unique strata (also works if 'breakdown' is NULL)
    ## create object of class "arpr"
    res <- constructRmpg(value=value, valueByStratum=valueByStratum, 
        years=ys, strata=rs, threshold=ts)
    # variance estimation (if requested)
    if(!is.null(var)) {
        res <- variance(inc, weights, years, breakdown, design, 
            indicator=res, alpha=alpha, na.rm=na.rm, type=var, ...)
    }
    ## return result
    return(res)
}

## workhorse
relativeGap <- function(x, weights = NULL, 
        sort = NULL, threshold, na.rm = FALSE) {
    ## initializations
    if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
    if(isTRUE(na.rm)){
        indices <- !is.na(x)
        x <- x[indices]
        if(!is.null(weights)) weights <- weights[indices]
        if(!is.null(sort)) sort <- sort[indices]
    } else if(any(is.na(x))) return(NA)
    if(length(x) == 0) return(NA)
    # preparations
    isPoor <- x < threshold  # individuals below threshold
    x <- x[isPoor]
    if(!is.null(weights)) weights <- weights[isPoor]
    if(!is.null(sort)) sort <- sort[isPoor]
    # calculations
    medianPoor <- incMedian(x, weights, sort)
    (threshold - medianPoor) * 100 / threshold
}
