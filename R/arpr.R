# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## 1. at-risk-of-poverty rate
arpr <- function(inc, weights = NULL, sort = NULL, years = NULL, 
        breakdown = NULL, design = NULL, data = NULL, p = 0.6, 
        var = NULL, alpha = 0.05, na.rm = FALSE, ...) {
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
    if(byYear) {  # ARPR by year
        ys <- sort(unique(years))
        ts <- arpt(inc, weights, sort, years, p=p, na.rm=na.rm)  # thresholds
        wr <- function(y, t, inc, weights, years, na.rm) {
            i <- years == y
            weightedRate(inc[i], weights[i], t, na.rm=na.rm)
        }
        value <- mapply(wr, y=ys, t=ts, 
            MoreArgs=list(inc=inc, weights=weights, years=years, na.rm=na.rm))
        names(value) <- ys  # use years as names
        if(byStratum) {
            wr1 <- function(i, inc, weights, years, ts, na.rm) {
                y <- years[i[1]]
                t <- ts[as.character(y)]
                weightedRate(inc[i], weights[i], t, na.rm=na.rm)
            }
            valueByStratum <- aggregate(1:n, list(year=years, stratum=breakdown), 
                wr1, inc=inc, weights=weights, years=years, ts=ts, na.rm=na.rm)
            names(valueByStratum)[3] <- "value"
        } else valueByStratum <- NULL
    } else {  # ARPR for only one year
        ys <- NULL
        ts <- arpt(inc, weights, sort, p=p, na.rm=na.rm)  # threshold
        value <- weightedRate(inc, weights, ts, na.rm=na.rm)
        if(byStratum) {
            wr2 <- function(i, inc, weights, ts, na.rm) {
                weightedRate(inc[i], weights[i], ts, na.rm=na.rm)
            }
            valueByStratum <- aggregate(1:n, list(stratum=breakdown), 
                wr2, inc=inc, weights=weights, ts=ts, na.rm=na.rm)
            names(valueByStratum)[2] <- "value"
        } else valueByStratum <- NULL
    }
    rs <- levels(breakdown)  # unique strata (also works if 'breakdown' is NULL)
    ## create object of class "arpr"
    res <- constructArpr(value=value, valueByStratum=valueByStratum, 
        years=ys, strata=rs, p=p, threshold=ts)
    # variance estimation (if requested)
    if(!is.null(var)) {
        res <- variance(inc, weights, years, breakdown, design, 
            indicator=res, alpha=alpha, na.rm=na.rm, type=var, ...)
    }
    # return results
    return(res)
}

## workhorse
weightedRate <- function(x, weights = NULL, threshold, na.rm = FALSE) {
    ## initializations
    if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
    if(isTRUE(na.rm)){
        indices <- !is.na(x)
        x <- x[indices]
        weights <- weights[indices]
    } else if(any(is.na(x))) return(NA)
    ## calculations
    sw <- sum(weights)  # estimate population total
#    sum(weights[x < threshold])/sw  # percentage of persons below threshold
    sum(weights[x < threshold])*100/sw  # percentage of persons below threshold
}
