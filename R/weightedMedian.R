# ------------------------------------------
# Authors: Andreas Alfons and Matthias Templ
#          Vienna University of Technology
# ------------------------------------------

# weighted median
weightedMedian <- function(x, weights = NULL, sorted = FALSE, na.rm = FALSE) {
    weightedQuantile(x, weights, probs=0.5, sorted=sorted, na.rm=na.rm)
}
