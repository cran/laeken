## gpg coefficient
gpg <- function(inc, gender = NULL, method = c("mean", "median"), 
		weights = NULL, sort = NULL, years = NULL, breakdown = NULL, 
		design = NULL, data = NULL, var = NULL, alpha = 0.05, 
		na.rm = FALSE, ...) {
    ## initializations
    if(is.null(gender)) stop("'gender' must be supplied")
    byYear <- !is.null(years)
    byStratum <- !is.null(breakdown)
    if(!is.null(data)) {
        inc <- data[, inc]
        gender <- data[, gender]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(sort)) sort <- data[, sort]
        if(byYear) years <- data[, years]
        if(byStratum) breakdown <- data[, breakdown]
        if(!is.null(var) && !is.null(design)) design <- data[, design]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
#    if(method != 'mean' && method != 'median') stop("'method' must be 'mean' or 'median'")
    method <- match.arg(method)
	if(!is.factor(gender)) stop("'gender' must be a factor.")
    if(length(levels(gender)) != 2) stop("'gender' must have exactly two levels: male = '1', female = '2'")
    if(!is.null(years)) {
    	if(!is.factor(years)) stop("'years' should be a factor")
	    nage <- length(levels(years))
    	if(n > 12) warning(paste("Too small sample sizes may occur by using ", n," age classes")) }
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
    # GPG by year (if requested)
    if(byYear) {
        ys <- sort(unique(years))  # unique years
        gc <- function(y, inc, weights, sort, years, na.rm) {
            i <- years == y
            genderGap(inc[i], gender[i], method, weights[i], sort[i], na.rm=na.rm)
        }
        value <- sapply(ys, gc, inc=inc, weights=weights, 
            sort=sort, years=years, na.rm=na.rm)
        names(value) <- ys  # use years as names
    } else {
        ys <- NULL
        value <- genderGap(inc, gender, method, weights, sort, na.rm=na.rm)
    }
    # GPG by stratum (if requested)
    if(byStratum) {
        gcR <- function(i, inc, weights, sort, na.rm) {
            genderGap(inc[i], gender[i], method, weights[i], sort[i], na.rm=na.rm)
        }
        valueByStratum <- aggregate(1:n, 
            if(byYear) list(year=years, stratum=breakdown) else list(stratum=breakdown), 
            gcR, inc=inc, weights=weights, sort=sort, na.rm=na.rm)
        names(valueByStratum)[ncol(valueByStratum)] <- "value"
        rs <- levels(breakdown)  # unique strata
    } else valueByStratum <- rs <- NULL
    ## create object of class "gpg"
    res <- constructGpg(value=value, 
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
genderGap <- function(x, gend, method = 'mean', weights = NULL, 
        sort = NULL, na.rm = FALSE) {
    
    if(is.null(gend)) stop("'gender' must be supplied")
    
    # initializations
    if(isTRUE(na.rm)){
        indices <- !is.na(x)
        x <- x[indices]
        gend <- gend[indices]
        if(!is.null(weights)) weights <- weights[indices]
        if(!is.null(sort)) sort <- sort[indices]
    } else if(any(is.na(x))) return(NA)
    
    male <- levels(gend)[1]
    female <- levels(gend)[2]
    
    
    if(is.null(weights)) weights <- rep.int(1, length(x))  # equal weights
    
    incgendmale <- x[gend=="male"]
    incgendmaleWeights <- weights[gend=="male"]
    incgendfemale <- x[gend=="female"]
    incgendfemaleWeights <- weights[gend=="female"]
    
    
    if(method == 'mean') {
        wM <- weighted.mean(x=incgendmale, w=incgendmaleWeights)
        wF <- weighted.mean(x=incgendfemale, w=incgendfemaleWeights)
        return((wM - wF) / wM)
        
    } else {
        wM <- weightedMedian(incgendmale, incgendmaleWeights)
        wF <- weightedMedian(incgendfemale, incgendfemaleWeights)
        return((wM - wF)/wM)
    }
}
