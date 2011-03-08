# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# get matrix of binary variables for calibration
calibVars <- function(x) UseMethod("calibVars")

calibVars.default <- function(x) {
    if(length(x) == 0) matrix(integer(), 0, 0)
    x <- as.factor(x)
    res <- sapply(levels(x), function(l) as.integer(x == l))
    rownames(res) <- names(x)  # set rownames from original vector
    res
}

calibVars.matrix <- function(x) calibVars(as.data.frame(x))

calibVars.data.frame <- function(x) {
    res <- lapply(x, calibVars)  # list of matrices for each variable
    res <- mapply(function(x, nam) {
            colnames(x) <- paste(nam, colnames(x), sep=".")
            x
        }, res, names(x), SIMPLIFY=FALSE)
    res <- do.call("cbind", res)  # combine matrices
    rownames(res) <- row.names(x)  # set rownames from original data.frame
    res
}

