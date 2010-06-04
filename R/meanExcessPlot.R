# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

meanExcessPlot <- function(x, probs, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(any(i <- is.na(x))) x <- x[!i]
    n <- length(x)
    if(n == 0) stop("no observed values")
    if(missing(probs)) {
        eps <- 1/n
        probs <- seq.int(from=eps, to=1-eps*sqrt(n), by=eps)
    }
    mu <- quantile(x, probs, type=1)
    meanExcess <- function(mu) {
        if(max(x) <= mu) stop("threshold too high")
        mean(x[x > mu] - mu)
    }
    me <- sapply(mu, meanExcess)
    ## create plot
    localPlot <- function(x, y, main = "Mean excess plot", 
            xlab = "Threshold", ylab = "Mean excess", ...) {
        plot(x, y, main=main, xlab=xlab, ylab=ylab, ...)
    }
    localPlot(mu, me, ...)
}
