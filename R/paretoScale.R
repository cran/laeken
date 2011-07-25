# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

paretoScale <- function(x, w = NULL, groups = NULL, 
        method = "VanKerm", center = c("mean", "median"), 
		probs = c(0.97, 0.98), na.rm = FALSE) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    useW <- !is.null(w)
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
        x <- x[unique]
        if(useW) w <- w[unique]
    }
#	method <- match.arg(method)  # only van Kerm's method currently implemented
	center <- match.arg(center)
	probs <- rep(probs, length.out=2)
    na.rm <- isTRUE(na.rm)
    # estimate threshold with van Kerm's method
    if(center == "mean") {
		mu <- weightedMean(x, w, na.rm=na.rm)
		q <- weightedQuantile(x, w, probs=probs, na.rm=na.rm)
	} else {
		q <- weightedQuantile(x, w, probs=c(0.5, probs), na.rm=na.rm)
		mu <- q[1]
		q <- q[-1]
	}
    x0 <- max(min(2.5*mu, q[2]), q[1])
    res <- list(x0=x0, k=length(which(x > x0)))
    class(res) <- "paretoScale"
    res
}


## print method for class "paretoScale"
print.paretoScale <- function(x, ...) {
    cat("Threshold: ")
    cat(x$x0, ...)
    cat("\nNumber of observations in the tail: ")
    cat(x$k, ...)
    cat("\n")
}
