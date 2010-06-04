# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

paretoQPlot <- function(x, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(any(i <- is.na(x))) x <- x[!i]
    x <- sort(x)
    n <- length(x)
    if(n == 0) stop("no observed values")
    y <- -log((n:1)/(n+1))
    ## create plot
    localPlot <- function(x, y, main = "Pareto quantile plot", 
            xlab = expression(paste(-log, " ", 
                    frac((n-i+1),(n+1)), " ,     i=1,...,n")), 
            ylab = expression(paste(log, " x")), log, xlog, ylog, ...) {
        plot(x, y, main=main, xlab=xlab, ylab=ylab, ...)
    }
    localPlot(y, log(x), ...)
}
