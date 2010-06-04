# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

variance <- function(inc, weights = NULL, years = NULL, breakdown = NULL, 
        design = NULL, data = NULL, indicator, alpha = 0.05, 
        na.rm = FALSE, type = "bootstrap", ...) {
    # initializations
    type <- match.arg(type)
    # call function corresponding to 'type'
    switch(type,
        bootstrap = bootVar(inc, weights, years, breakdown, design, 
            data, indicator, alpha=alpha, na.rm=na.rm, ...))
}
