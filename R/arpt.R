# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# at-risk-of poverty threshold
arpt <- function(inc, weights = NULL, sort = NULL, 
        years = NULL, data = NULL, p = 0.6, na.rm = FALSE) {
    # check 'p' (other arguments are checked in 'incMedian')
    if(!is.numeric(p) || length(p) == 0 || p[1] < 0 || p[1] > 1) {
        stop("'p' must be a numeric value in [0,1]")
    } else p <- p[1]
    # compute at-risk-of-poverty threshold
    p * incMedian(inc, weights, sort, years, data, na.rm=na.rm)
}
