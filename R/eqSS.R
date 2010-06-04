# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# TODO: error handling

# equivalised household size
eqSS <- function(hid, age, year = NULL, data = NULL) {
    ## initializations
    if(is.null(data)) {
        data <- data.frame(hid=hid)
        hid <- "hid"
        if(!is.null(year)) {
            data <- cbind(year=year, data)
            year <- "year"
        }
    } else {
        age <- data[, age]
        data <- data[, c(year, hid), drop=FALSE]
    }
    ## calculations
    i <- if(is.null(year)) 2 else 3
    tmp <- as.data.frame(table(data))  # number of household members
    hm14p <- as.data.frame(table(data[age >= 14,]))[, i]  # at least 14 years
    hm13m <- tmp[, i] - hm14p  # younger than 14
    tmp[, i] <- 1 + 0.5*(hm14p-1) + 0.3*hm13m  # eqSS for househoulds
    names(tmp) <- c(year, hid, ".eqSS")
    data <- cbind(data, .ID=1:nrow(data))  # add ID to original data
    data <- merge(data, tmp, sort=FALSE)  # merge with original data set
    ## order according to original data and extract eqSS
    data$.eqSS[order(data$.ID)]
}
