# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## FIXME: do not use 'p' as argument name for function passed to 'boot'

## generic function
bootVar <- function(inc, weights = NULL, years = NULL, 
		breakdown = NULL, design = NULL, data = NULL, indicator, 
		R = 100, bootType = c("calibrate", "naive"), X, 
		totals = NULL, ciType = c("perc", "norm", "basic"),
		# type "stud" and "bca" are currently not allowed
		alpha = 0.05, seed = NULL, na.rm = FALSE, gender = NULL, method = 'mean', ...) {
	UseMethod("bootVar", indicator)
}


## class "indicator"
bootVar.indicator <- function(inc, weights = NULL, years = NULL, 
		breakdown = NULL, design = NULL, data = NULL, indicator, 
		R = 100, bootType = c("calibrate", "naive"), X, 
		totals = NULL, ciType = c("perc", "norm", "basic"),
		# type "stud" and "bca" are currently not allowed
		alpha = 0.05, seed = NULL, na.rm = FALSE, gender = NULL, method = 'mean', ...) {
	## initializations
	# check whether weights have been supplied
	haveWeights <- !is.null(weights)    
	haveGender <- !is.null(gender) 
	# check whether indicator is broken down by year
	# if so, check whether years have been supplied
	ys <- indicator$years
	byYear <- !is.null(ys)
	if(byYear && is.null(years)) stop("'years' must be supplied")
	# check whether indicator is broken down by stratum
	# if so, check whether breakdown has been supplied
	rs <- indicator$strata
	byStratum <- !is.null(rs)
	if(byStratum && is.null(breakdown)) stop("'breakdown' must be supplied")
	haveDesign <- !is.null(design)
	# if a data.frame has been supplied, extract the respective vectors
	if(!is.null(data)) {
		inc <- data[, inc]
		if(!is.null(weights)) weights <- data[, weights]
		if(!is.null(gender)) gender <- data[, gender]
		if(byYear) years <- data[, years]
		if(byStratum) breakdown <- data[, breakdown]
		if(haveDesign) design <- data[, design]
	}
	# check whether the vectors have the correct type
	if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
	n <- length(inc)
	if(haveWeights && !is.numeric(weights)) {
		stop("'weights' must be a numeric vector")
	}
#	if(haveGender && !is.numeric(gender)) {
#		stop("'gender' must be a numeric vector")
#	}
	if(byYear && !is.numeric(years)) {
		stop("'years' must be a numeric vector")
	}
	if(byStratum && !is.vector(breakdown) && !is.factor(breakdown)) {
		stop("'breakdown' must be a vector or factor")
	}
	if(haveDesign && !is.integer(design) && !is.factor(design)) {
		stop("'design' must be an integer vector or factor")
	}
	if(is.null(data)) {  # check vector lengths
		if(haveWeights && length(weights) != n) {
			stop("'weights' must have length ", n)
		}
		if(byYear && length(years) != n) {
			stop("'years' must have length ", n)
		}
		if(byStratum && length(breakdown) != n) {
			stop("'breakdown' must have length ", n)
		}
		if(haveDesign && length(design) != n) {
			stop("'design' must have length ", n)
		}
	}
	if(!haveDesign) design <- rep.int(1, n)
	# check other input
	if(!is.numeric(R) || length(R) == 0) stop("'R' must be numeric")
	else R <- as.integer(R[1])
	if(!is.numeric(alpha) || length(alpha) == 0) stop("'alpha' must be numeric")
	else alpha <- alpha[1]
	bootType <- match.arg(bootType)
	calibrate <- haveWeights && bootType == "calibrate"
	if(calibrate) {
		X <- as.matrix(X)
#        if(!is.numeric(X)) stop("'X' must be a numeric matrix")
		if(nrow(X) != n) stop("'X' must have ", n, " rows")
		if(is.null(totals)) {
			# compute totals from original data with Horvitz-Thompson estimator
			if(byYear) {
				totals <- lapply(ys, 
						function(y) {
							# extract current year from calibration variables and 
							# weights
							i <- years == y
							X <- X[i, , drop=FALSE]
							weights <- weights[i]
							# compute totals for current year
							apply(X, 2, function(i) sum(i*weights))
						})
				totals <- do.call(rbind, totals)  # form matrix of totals
				rownames(totals) <- ys  # use years as rownames for totals
			} else totals <- apply(X, 2, function(i) sum(i*weights))
		} else if(byYear) totals <- as.matrix(totals)
		if(!is.numeric(totals)) stop("'totals' must be of type numeric")
	} else {
		X <- NULL
		totals <- NULL
	}
	ciType <- match.arg(ciType)
	
	## preparations
	data <- data.frame(inc=inc)
	data$gender <- gender
	data$weight <- weights
	data$year <- years
	data$stratum <- breakdown
	data$method <- method
	if(inherits(indicator, "arpr")) {
		p <- indicator$p  # percentage of median used for threshold
	} else p <- NULL
	if(!is.null(seed)) set.seed(seed)  # set seed of random number generator
	seed <- .Random.seed  # seed is later on added to the object to be returned
	
	## calculations
	# get basic function for bootstrap replications with definition:
	# function(x, i, p, X, totals, rs, na.rm)
	fun <- getFun(indicator, byStratum)
	bootFun <- getBootFun(calibrate, fun)
	if(byYear) {
		# ---------- breakdown by year ----------
		# get more complex function for additional with definition
		# function(y, x, R, p, aux, totals, rs, alpha, ciType, na.rm, ...)
		funByYear <- getFunByYear(byStratum, calibrate, bootFun)
		if(byStratum) {
			# ---------- breakdown by stratum ----------
			tmp <- lapply(ys, funByYear, data, R, design, p, X, 
					totals, ys, rs, alpha, ciType, na.rm, ...)
			var <- do.call(c, lapply(tmp, function(x) x[[1]]))
			names(var) <- ys
			varByStratum <- do.call(rbind, lapply(tmp, function(x) x[[2]]))
			ci <- do.call(rbind, lapply(tmp, function(x) x[[3]]))
			rownames(ci) <- ys
			ciByStratum <- do.call(rbind, lapply(tmp, function(x) x[[4]]))
			# order 'varByStratum' and 'ciByStratum' according to 'valueByStratum'
			tmp <- indicator$valueByStratum[, 1:2]
			tmp <- data.frame(order=1:nrow(tmp), tmp)
			varByStratum <- merge(varByStratum, tmp, all=TRUE, sort=FALSE)
			varByStratum <- varByStratum[order(varByStratum$order), -4]
			ciByStratum <- merge(ciByStratum, tmp, all=TRUE, sort=FALSE)
			ciByStratum <- ciByStratum[order(ciByStratum$order), -5]
		} else {
			# ---------- no breakdown by stratum ----------
			tmp <- sapply(ys, funByYear, data, R, design, p, X, 
					totals, ys, rs, alpha, ciType, na.rm, ...)
			colnames(tmp) <- ys
			var <- tmp[1,]
			ci <- t(tmp[2:3,])
		}
	} else {
		# ---------- no breakdown by year ----------
		b <- boot(data, bootFun, R, strata=design, p=p, aux=X, 
				totals=totals, rs=rs, na.rm=na.rm, ...)
		if(byStratum) {
			# ---------- breakdown by stratum ----------
			var <- apply(b$t, 2, var)
			varByStratum <- data.frame(stratum=rs, var=var[-1])
			var <- var[1]
			ci <- sapply(1:((length(rs) + 1)), 
					function(i) {
						ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
						switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
								basic=ci$basic[4:5], stud=ci$student[4:5], 
								bca=ci$bca[4:5])
					})
			rownames(ci) <- c("lower", "upper")
			ciByStratum <- data.frame(stratum=rs, t(ci[, -1]))
			ci <- ci[, 1]
		} else {
			# ---------- no breakdown by stratum ----------
			var <- var(b$t[, 1])
			ci <- boot.ci(b, conf=1-alpha, type=ciType)
			ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
					basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
			names(ci) <- c("lower", "upper")
		}
	}
	
	## modify and return object
	indicator$varMethod <- "bootstrap"
	indicator$var <- var
	indicator$ci <- ci
	if(byStratum) {
		indicator$varByStratum <- varByStratum
		indicator$ciByStratum <- ciByStratum
	}
	indicator$alpha <- alpha
	indicator$seed <- seed
	return(indicator)
	#return(data$method)
}


## utility functions: return functions to be used in the bootstrap replications


# basic function for breakdown by stratum
getFun <- function(indicator, byStratum) UseMethod("getFun")

getFun.arpr <- function(indicator, byStratum) {
	if(byStratum) {
		function(x, p, rs, na.rm) {
			threshold <- p * weightedMedian(x$inc, x$weight)
			value <- weightedRate(x$inc, x$weight, threshold, na.rm=na.rm)
			valueByStratum <- sapply(rs, function(r, x, t) {
						i <- x$stratum == r
						weightedRate(x$inc[i], x$weight[i], t, na.rm=na.rm)
					}, x=x, t=threshold)
			c(value, valueByStratum)
		}
	} else {
		function(x, p, rs, na.rm) {
			threshold <- p * weightedMedian(x$inc, x$weight)
			weightedRate(x$inc, x$weight, threshold, na.rm=na.rm)
		}
	}
}

# the argument 'p' is not necessary here, but is used so 
# that we have a unified function call for all indicators
getFun.qsr <- function(indicator, byStratum) {
	if(byStratum) {
		function(x, p, rs, na.rm) {
			value <- quintileRatio(x$inc, x$weight, na.rm=na.rm)
			valueByStratum <- sapply(rs, function(r, x, t) {
						i <- x$stratum == r
						quintileRatio(x$inc[i], x$weight[i], na.rm=na.rm)
					}, x=x)
			c(value, valueByStratum)
		}
	} else {
		function(x, p, rs, na.rm) {
			quintileRatio(x$inc, x$weight, na.rm=na.rm)
		}
	}
}

getFun.rmpg <- function(indicator, byStratum) {
	if(byStratum) {
		function(x, p, rs, na.rm) {
			threshold <- 0.6 * weightedMedian(x$inc, x$weight)
			value <- relativeGap(x$inc, x$weight, 
					threshold=threshold, na.rm=na.rm)
			valueByStratum <- sapply(rs, function(r, x, t) {
						i <- x$stratum == r
						relativeGap(x$inc[i], x$weight[i], threshold=t, na.rm=na.rm)
					}, x=x, t=threshold)
			c(value, valueByStratum)
		}
	} else {
		function(x, p, rs, na.rm) {
			threshold <- 0.6 * weightedMedian(x$inc, x$weight)
			relativeGap(x$inc, x$weight, threshold=threshold, na.rm=na.rm)
		}
	}
}

# the argument 'p' is not necessary here, but is used so 
# that we have a unified function call for all indicators
getFun.gini <- function(indicator, byStratum) {
	if(byStratum) {
		function(x, p, rs, na.rm) {
			value <- giniCoeff(x$inc, x$weight, na.rm=na.rm)
			valueByStratum <- sapply(rs, function(r, x, t) {
						i <- x$stratum == r
						giniCoeff(x$inc[i], x$weight[i], na.rm=na.rm)
					}, x=x)
			c(value, valueByStratum)
		}
	} else {
		function(x, p, rs, na.rm) {
			giniCoeff(x$inc, x$weight, na.rm=na.rm)
		}
	}
}

# the argument 'p' is not necessary here, but is used so 
# that we have a unified function call for all indicators
getFun.gpg <- function(indicator, byStratum) {
	if(byStratum) {
		function(x, p, rs, na.rm) {
			value <- gpgCoeff(x$inc, x$gender, x$method[1], x$weight, na.rm=na.rm)
			valueByStratum <- sapply(rs, function(r, x, t) {
						i <- x$stratum == r
						gpgCoeff(x$inc[i], x$gender[i], x$method[1], x$weight[i], na.rm=na.rm)
					}, x=x)
			c(value, valueByStratum)
		}
	} else {
		function(x, p, rs, na.rm) {
			gpgCoeff(x$inc, x$gender, x$method[1], x$weight, na.rm=na.rm)
		}
	}
}


# function that incorporates resampling and (if requested) calibration
getBootFun <- function(calibrate, fun) {
	if(calibrate) {
		function(x, i, p, aux, totals, rs, na.rm, ...) {
			x <- x[i, , drop=FALSE]
			aux <- aux[i, , drop=FALSE]
			g <- calibWeights(aux, x$weight, totals, ...)
			x$weight <- g * x$weight
			fun(x, p, rs, na.rm)
		}
	} else {
		function(x, i, p, aux, totals, rs, na.rm, ...) {
			x <- x[i, , drop=FALSE]
			fun(x, p, rs, na.rm)
		}
	}
}


# more complex function for additional breakdown by year
getFunByYear <- function(byStratum, calibrate, fun) {
	if(byStratum) {
		if(calibrate) {
			# ---------- breakdown by stratum, calibration ----------
			function(y, x, R, design, p, aux, totals, ys, rs, alpha, ciType, na.rm, ...) {
				i <- x$year == y
				x <- x[i, , drop=FALSE]
				aux <- aux[i, , drop=FALSE]
				design <- design[i]
				totals <- totals[ys == y,]
				b <- boot(x, fun, R, strata=design, p=p, aux=aux, 
						totals=totals, rs=rs, na.rm=na.rm, ...)
				var <- apply(b$t, 2, var)
				varByStratum <- data.frame(year=y, stratum=rs, var=var[-1])
				var <- var[1]
				ci <- sapply(1:((length(rs) + 1)), 
						function(i) {
							ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
							switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
									basic=ci$basic[4:5], stud=ci$student[4:5], 
									bca=ci$bca[4:5])
						})
				rownames(ci) <- c("lower", "upper")
				ciByStratum <- data.frame(year=y, stratum=rs, t(ci[, -1]))
				ci <- ci[, 1]
				list(var, varByStratum, ci, ciByStratum)
			}
		} else {
			# ---------- breakdown by stratum, no calibration ----------
			function(y, x, R, design, p, aux, totals, ys, rs, alpha, ciType, na.rm, ...) {
				i <- x$year == y
				x <- x[i, , drop=FALSE]
				design <- design[i]
				b <- boot(x, fun, R, strata=design, p=p, aux=aux, 
						totals=totals, rs=rs, na.rm=na.rm, ...)
				var <- apply(b$t, 2, var)
				varByStratum <- data.frame(year=y, stratum=rs, var=var[-1])
				var <- var[1]
				ci <- sapply(1:((length(rs) + 1)), 
						function(i) {
							ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
							switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
									basic=ci$basic[4:5], stud=ci$student[4:5], 
									bca=ci$bca[4:5])
						})
				rownames(ci) <- c("lower", "upper")
				ciByStratum <- data.frame(year=y, stratum=rs, t(ci[, -1]))
				ci <- ci[, 1]
				list(var, varByStratum, ci, ciByStratum)
			}
		}
	} else {
		if(calibrate) {
			# ---------- no breakdown by stratum, calibration ----------
			function(y, x, R, design, p, aux, totals, ys, rs, alpha, ciType, na.rm, ...) {
				i <- x$year == y
				x <- x[i, , drop=FALSE]
				aux <- aux[i, , drop=FALSE]
				design <- design[i]
				totals <- totals[ys == y,]
				b <- boot(x, fun, R, strata=design, p=p, aux=aux, 
						totals=totals, rs=rs, na.rm=na.rm, ...)
				var <- var(b$t[, 1])
				ci <- boot.ci(b, conf=1-alpha, type=ciType)
				ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
						basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
				names(ci) <- c("lower", "upper")
				c(var, ci)
			}
		} else {
			# ---------- no breakdown by stratum, no calibration ----------
			function(y, x, R, design, p, aux, totals, ys, rs, alpha, ciType, na.rm, ...) {
				i <- x$year == y
				x <- x[i, , drop=FALSE]
				design <- design[i]
				b <- boot(x, fun, R, strata=design, p=p, aux=aux, 
						totals=totals, rs=rs, na.rm=na.rm, ...)
				var <- var(b$t[, 1])
				ci <- boot.ci(b, conf=1-alpha, type=ciType)
				ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3], 
						basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
				names(ci) <- c("lower", "upper")
				c(var, ci)
			}
		}
	}
}
