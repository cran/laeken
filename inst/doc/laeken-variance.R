### R code from vignette source 'laeken-variance.Rnw'

###################################################
### code chunk number 1: laeken-variance.Rnw:52-53
###################################################
options(prompt="R> ")


###################################################
### code chunk number 2: laeken-variance.Rnw:93-95 (eval = FALSE)
###################################################
## vignette("laeken-standard")
## vignette("laeken-pareto")


###################################################
### code chunk number 3: laeken-variance.Rnw:111-113
###################################################
library("laeken")
data("eusilc")


###################################################
### code chunk number 4: laeken-variance.Rnw:140-141
###################################################
args(variance)


###################################################
### code chunk number 5: laeken-variance.Rnw:226-229
###################################################
a <- arpr("eqIncome", weights = "rb050", data = eusilc)
variance("eqIncome", weights = "rb050", design = "db040",
    data = eusilc, indicator = a, bootType = "naive", seed = 123)


###################################################
### code chunk number 6: laeken-variance.Rnw:237-240
###################################################
b <- arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)
variance("eqIncome", weights = "rb050", breakdown = "db040", design = "db040",
    data = eusilc, indicator = b, bootType = "naive", seed = 123)


###################################################
### code chunk number 7: laeken-variance.Rnw:309-312
###################################################
variance("eqIncome", weights = "rb050", design = "db040",
    data = eusilc, indicator = a, X = calibVars(eusilc$db040),
    seed = 123)


