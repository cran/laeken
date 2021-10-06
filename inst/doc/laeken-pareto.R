### R code from vignette source 'laeken-pareto.Rnw'

###################################################
### code chunk number 1: laeken-pareto.Rnw:74-75
###################################################
options(prompt="R> ")


###################################################
### code chunk number 2: laeken-pareto.Rnw:149-151 (eval = FALSE)
###################################################
## vignette("laeken-standard")
## vignette("laeken-variance")


###################################################
### code chunk number 3: laeken-pareto.Rnw:165-167
###################################################
library("laeken")
data("eusilc")


###################################################
### code chunk number 4: laeken-pareto.Rnw:235-236
###################################################
qsr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 5: laeken-pareto.Rnw:260-261
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 6: laeken-pareto.Rnw:293-298
###################################################
x <- seq(1, 6, length.out=1000)
dpareto <- function(x, x0 = 1, theta = 1) theta*x0^theta / x^(theta+1)
y1 <- dpareto(x, theta=1)
y2 <- dpareto(x, theta=2)
y3 <- dpareto(x, theta=3)


###################################################
### code chunk number 7: laeken-pareto.Rnw:303-313
###################################################
par(mar = c(4, 4, 0.5, 0.5) + 0.1)
plot(x, y3, type = "l", lty = 3, ylab = "f(x)", xlim = c(0.75, 6),
    panel.first = {
        abline(h = 0, col = grey(0.75))
        abline(v = 1, col = grey(0.75))
    })
lines(x, y2, lty = 2)
lines(x, y1, lty = 1)
leg <- expression(paste(theta, " = 1"), paste(theta, " = 2"), paste(theta, " = 3"))
legend("topright", legend = leg, lty = 1:3)


###################################################
### code chunk number 8: laeken-pareto.Rnw:355-357
###################################################
hID <- eusilc$db030[which.max(eusilc$eqIncome)]
eusilc[eusilc$db030 == hID, "eqIncome"] <- 10000000


###################################################
### code chunk number 9: laeken-pareto.Rnw:366-367
###################################################
eusilcH <- eusilc[!duplicated(eusilc$db030), c("eqIncome", "db090")]


###################################################
### code chunk number 10: laeken-pareto.Rnw:424-426
###################################################
ts <- paretoScale(eusilcH$eqIncome, w = eusilcH$db090)
ts


###################################################
### code chunk number 11: laeken-pareto.Rnw:491-492
###################################################
paretoQPlot(eusilcH$eqIncome, w = eusilcH$db090)


###################################################
### code chunk number 12: laeken-pareto.Rnw:539-540
###################################################
meanExcessPlot(eusilcH$eqIncome, w = eusilcH$db090)


###################################################
### code chunk number 13: laeken-pareto.Rnw:592-594
###################################################
thetaHill(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaHill(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 14: laeken-pareto.Rnw:671-673
###################################################
thetaWML(eusilcH$eqIncome, k = ts$k)
thetaWML(eusilcH$eqIncome, x0 = ts$x0)


###################################################
### code chunk number 15: laeken-pareto.Rnw:718-720
###################################################
thetaISE(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaISE(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 16: laeken-pareto.Rnw:760-762
###################################################
thetaPDC(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaPDC(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 17: laeken-pareto.Rnw:800-801
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 18: laeken-pareto.Rnw:812-814
###################################################
fit <- paretoTail(eusilc$eqIncome, k = ts$k,
    w = eusilc$db090, groups = eusilc$db030)


###################################################
### code chunk number 19: laeken-pareto.Rnw:824-826
###################################################
w <- reweightOut(fit, calibVars(eusilc$db040))
gini(eusilc$eqIncome, w)


###################################################
### code chunk number 20: laeken-pareto.Rnw:834-837
###################################################
set.seed(1234)
eqIncome <- replaceOut(fit)
gini(eqIncome, weights = eusilc$rb050)


###################################################
### code chunk number 21: laeken-pareto.Rnw:842-845
###################################################
set.seed(1234)
eqIncome <- replaceTail(fit)
gini(eqIncome, weights = eusilc$rb050)


