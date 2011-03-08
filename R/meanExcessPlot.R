# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

meanExcessPlot <- function(x, w = NULL, probs = NULL, 
        interactive = TRUE, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(is.null(w)) {  # no weights
        if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
    } else {  # weights are supplied
        if(!is.numeric(w) || length(w) != length(x)) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
        if(any(i <- is.na(x))) {  # remove missing values
            x <- x[!i]
            w <- w[!i]
        }
    }
    n <- length(x)
    if(n == 0) stop("no observed values")
    ## use observed values or quantiles as thresholds
    if(is.null(probs)) {
#        eps <- 1/n
#        probs <- seq.int(from=eps, to=1-eps*sqrt(n), by=eps)
        mu <- unname(sort(x)[1:(n-sqrt(n))])
    } else {
        if(is.null(w)) {  # no weights
            mu <- quantile(x, probs, names=FALSE, type=1)  # compute quantiles
        } else {  # weights are supplied
            mu <- weightedQuantile(x, probs)  # compute weighted quantiles
        }
        if(max(mu) >= max(x)) stop("largest threshold too high")
    }
    ## compute mean excesses for the different thresholds
    # this could be done much faster with C (incremental computation)
#    meanExcess <- function(mu) {
#        if(max(x) <= mu) stop("threshold too high")
#        mean(x[x > mu] - mu)
#    }
    if(is.null(w)) meanExcess <- function(mu) mean(x[x > mu] - mu)
    else {
        meanExcess <- function(mu) {
            i <- x > mu
            weighted.mean(x[i] - mu, w[i])
        }
    }
    me <- sapply(mu, meanExcess)
    ## create plot
    localPlot <- function(x, y, main = "Mean excess plot", 
            xlab = "Threshold", ylab = "Mean excess", ...) {
        plot(x, y, main=main, xlab=xlab, ylab=ylab, ...)
    }
    localPlot(mu, me, ...)
    ## interactive identification of threshold
    res <- NULL
    if(isTRUE(interactive)) {
        index <- identify(mu, me, n=1, plot=FALSE)
        i <- 1
        while(!identical(index, integer())) {
            last <- index
            x0 <- mu[index]
            res <- list(x0=x0, k=length(which(x > x0)))
            class(res) <- "paretoScale"
            if(i > 1) cat("\n")
            print(res)
            index <- identify(mu, me, n=1, plot=FALSE)
            i <- i + 1
        }
        ## plot horizontal and vertical lines for selected threshold
        if(!is.null(res)) {
            abline(h=me[last], lty=3)
            abline(v=x0, lty=3)
        }
    }
    ## return result invisibly
    invisible(res)
}
