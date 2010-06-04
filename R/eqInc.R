# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# TODO: error handling
# TODO: account for inflation-adjustment

# equivalised disposable income
eqInc <- function(hid, hplus, hminus, pplus, pminus, 
        eqSS, year = NULL, data = NULL) {
    ## initializations
    if(is.null(data)) {
        data <- data.frame(hid=hid)
        hid <- "hid"
        if(!is.null(year)) {
            data <- cbind(year=year, data)
            year <- "year"
        }
        npplus <- names(pplus)
        npminus <- names(pminus)
    } else {
        hplus <- data[, hplus]
        hminus <- data[, hminus]
        npplus <- pplus
        pplus <- data[, npplus]
        npminus <- pminus
        pminus <- data[, npminus]
        eqSS <- data[, eqSS]
        data <- data[, c(year, hid), drop=FALSE]
    }
    ## calculations
    hy020h <- rowSums(hplus, na.rm=TRUE) - rowSums(hminus, na.rm=TRUE)
    tmp <- aggregate(data.frame(pplus,pminus), data, sum, na.rm=TRUE)
    hy020p <- rowSums(tmp[,npplus], na.rm=TRUE) - 
        rowSums(tmp[,npminus], na.rm=TRUE)
    if(is.null(year)) {
        names(hy020p) <- tmp[, hid]
        hy020p <- unname(hy020p[as.character(data[, hid])])
    } else {
        tmp <- cbind(tmp[, c(year, hid), drop=FALSE], .hy020p=hy020p)
        data <- cbind(data, .ID=1:nrow(data))  # add ID to original data
        data <- merge(data, tmp, sort=FALSE)  # merge with original data set
        ## order according to original data and extract hy020p
        hy020p <- data$.hy020p[order(data$.ID)]
    }
    ## return result
    (hy020h + hy020p) / eqSS
}
