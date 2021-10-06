### R code from vignette source 'laeken-intro.Rnw'

###################################################
### code chunk number 1: laeken-intro.Rnw:107-109
###################################################
options(prompt = "R> ", continue = "+  ", width = 72, useFancyQuotes = FALSE)
library("laeken")


###################################################
### code chunk number 2: laeken-intro.Rnw:164-165 (eval = FALSE)
###################################################
## vignette(package="laeken")


###################################################
### code chunk number 3: laeken-intro.Rnw:244-246
###################################################
data("eusilc")
head(eusilc[, 1:10], 3)


###################################################
### code chunk number 4: laeken-intro.Rnw:269-271
###################################################
data("ses")
head(ses[, 1:7], 3)


###################################################
### code chunk number 5: laeken-intro.Rnw:392-393
###################################################
arpr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 6: laeken-intro.Rnw:408-409
###################################################
arpr("eqIncome", weights = "rb050", p = c(0.4, 0.5, 0.7), data = eusilc)


###################################################
### code chunk number 7: laeken-intro.Rnw:431-432
###################################################
qsr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 8: laeken-intro.Rnw:462-463
###################################################
rmpg("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 9: laeken-intro.Rnw:483-484
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 10: laeken-intro.Rnw:526-527
###################################################
gpg("earningsHour", gender = "sex", weigths = "weights", data = ses)


###################################################
### code chunk number 11: laeken-intro.Rnw:550-552
###################################################
gpg("earningsHour", gender = "sex", weigths = "weights", data = ses,
    method = "median")


###################################################
### code chunk number 12: laeken-intro.Rnw:589-590
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 13: laeken-intro.Rnw:593-594
###################################################
gini(eusilc$eqIncome, weights = eusilc$rb050)


###################################################
### code chunk number 14: laeken-intro.Rnw:670-672
###################################################
a <- arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)
a


###################################################
### code chunk number 15: laeken-intro.Rnw:686-687
###################################################
subset(a, strata = c("Lower Austria", "Vienna"))


###################################################
### code chunk number 16: laeken-intro.Rnw:755-758
###################################################
hID <- eusilc$db030[which.max(eusilc$eqIncome)]
eqIncomeOut <- eusilc$eqIncome
eqIncomeOut[eusilc$db030 == hID] <- 10000000


###################################################
### code chunk number 17: laeken-intro.Rnw:765-767
###################################################
keep <- !duplicated(eusilc$db030)
eusilcH <- data.frame(eqIncome=eqIncomeOut, db090=eusilc$db090)[keep,]


###################################################
### code chunk number 18: laeken-intro.Rnw:796-797
###################################################
paretoQPlot(eusilcH$eqIncome, w = eusilcH$db090)


###################################################
### code chunk number 19: laeken-intro.Rnw:852-854
###################################################
ts <- paretoScale(eusilcH$eqIncome, w = eusilcH$db090)
ts


###################################################
### code chunk number 20: laeken-intro.Rnw:919-921
###################################################
thetaISE(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaISE(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 21: laeken-intro.Rnw:953-955
###################################################
thetaPDC(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaPDC(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 22: laeken-intro.Rnw:1009-1011
###################################################
fit <- paretoTail(eqIncomeOut, k = ts$k, w = eusilc$db090,
                  groups = eusilc$db030)


###################################################
### code chunk number 23: laeken-intro.Rnw:1023-1024
###################################################
plot(fit)


###################################################
### code chunk number 24: laeken-intro.Rnw:1050-1052
###################################################
w <- reweightOut(fit, calibVars(eusilc$db040))
gini(eqIncomeOut, w)


###################################################
### code chunk number 25: laeken-intro.Rnw:1060-1063
###################################################
set.seed(123)
eqIncomeRN <- replaceOut(fit)
gini(eqIncomeRN, weights = eusilc$rb050)


###################################################
### code chunk number 26: laeken-intro.Rnw:1069-1071
###################################################
eqIncomeSN <- shrinkOut(fit)
gini(eqIncomeSN, weights = eusilc$rb050)


###################################################
### code chunk number 27: laeken-intro.Rnw:1078-1079
###################################################
gini(eqIncomeOut, weights = eusilc$rb050)


###################################################
### code chunk number 28: laeken-intro.Rnw:1152-1154
###################################################
arpr("eqIncome", weights = "rb050", design = "db040", cluster = "db030",
     data = eusilc, var = "bootstrap", bootType = "naive", seed = 1234)


###################################################
### code chunk number 29: laeken-intro.Rnw:1202-1205
###################################################
aux <- cbind(calibVars(eusilc$db040), calibVars(eusilc$rb090))
arpr("eqIncome", weights = "rb050", design = "db040", cluster = "db030",
     data = eusilc, var = "bootstrap", X = aux, seed = 1234)


