# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

paretoQPlot <- function(x, w = NULL, xlab = NULL, ylab = NULL, 
        interactive = TRUE, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(is.null(w)) {  # no weights
        if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
        x <- sort(x)  # sort values
    } else {  # weights are supplied
        if(!is.numeric(w) || length(w) != length(x)) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
        if(any(i <- is.na(x))) {  # remove missing values
            x <- x[!i]
            w <- w[!i]
        }
        # sort values and weights
        order <- order(x)
        x <- x[order]
        w <- w[order]
    }
    n <- length(x)
    if(n == 0) stop("no observed values")
    ## computation of theoretical quantiles
    if(is.null(w)) {
        y <- -log((n:1)/(n+1))
#        xlab <- expression(paste(-log, " ", 
#                frac((n+1-i),(n+1)), " ,     i=1,...,n"))
    } else {
        cw <- cumsum(w)
        y <- -log(1 - cw/(cw[n]*(n+1)/n))
#        xlab <- ""
    }
    ## create plot
#    ylab <- expression(paste(log, " x"))
    if(is.null(xlab)) xlab <- "Theoretical quantiles"
    if(is.null(ylab)) ylab <- ""
    localPlot <- function(x, y, main = "Pareto quantile plot", 
            xlab = NULL, ylab = NULL, log, xlog, ylog, ...) {
        plot(x, y, main=main, xlab=xlab, ylab=ylab, log="y", ...)
    }
#    localPlot(y, log(x), xlab=xlab, ylab=ylab, ...)
    localPlot(y, x, xlab=xlab, ylab=ylab, ...)
    ## interactive identification of threshold
    res <- NULL
    if(isTRUE(interactive)) {
#        index <- identify(y, log(x), n=1, plot=FALSE)
        index <- identify(y, x, n=1, plot=FALSE)
        i <- 1
        while(!identical(index, integer())) {
            last <- index
            x0 <- unname(x[index])
            res <- list(x0=x0, k=length(which(x > x0)))
            class(res) <- "paretoScale"
            if(i > 1) cat("\n")
            print(res)
#            index <- identify(y, log(x), n=1, plot=FALSE)
            index <- identify(y, x, n=1, plot=FALSE)
            i <- i + 1
        }
        ## plot horizontal and vertical lines for selected threshold
        if(!is.null(res)) {
#            abline(h=log(x0), lty=3)
            abline(h=x0, lty=3)
            abline(v=y[last], lty=3)
        }
    }
    ## return result invisibly
    invisible(res)
}
