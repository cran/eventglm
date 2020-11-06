## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 7, fig.height = 4.75
)

## ----setup--------------------------------------------------------------------
library(survival)
library(eventglm)

## ---- eval = FALSE------------------------------------------------------------
#  ?eventglm::colon

## -----------------------------------------------------------------------------
sfit <- survfit(Surv(time, status) ~ rx, data = colon)
plot(sfit, col = c("black", "slateblue", "salmon"), 
     xlab = "days since registration", ylab = "survival")
legend("bottomleft", fill = c("black", "slateblue", "salmon"), 
       legend = names(sfit$strata))

## -----------------------------------------------------------------------------
plot(sfit[1], conf.int = FALSE, xlab = "days since registration", ylab = "survival")

seg0 <- summary(sfit[1], times = sfit[1]$time[sfit[1]$time <= 2500])
rect(c(0, seg0$time), 0, c(seg0$time, 2500), c(seg0$surv), 
     border = NA, col = "grey80")
lines(sfit[1], conf.int = FALSE)
abline(v = 2500, lty = 2)
points(x = 2500, y = summary(sfit[1], times = 2500)$surv)


## -----------------------------------------------------------------------------
colon.sfit <- summary(sfit, times = 2500, rmean = 2500)
colon.sfit

## -----------------------------------------------------------------------------
colon.cifit <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon)
summary(colon.cifit)
se.ci <- sqrt(diag(vcov(colon.cifit, type = "robust")))
b.ci <- colon.cifit$coefficients
conf.ci <- confint(colon.cifit)

## -----------------------------------------------------------------------------
cbind(eventglm = b.ci, 
      survfit = c(1 - colon.sfit$surv[1], 
  (1 - colon.sfit$surv[2:3]) - 
    (1 - rep(colon.sfit$surv[1], 2))))

## -----------------------------------------------------------------------------
colon.rr <- cumincglm(Surv(time, status) ~ rx, time = 2500, 
                      data = colon, link = "log")
br.ci <- colon.rr$coefficients
confr.ci <- confint(colon.rr)

## -----------------------------------------------------------------------------
colon.rmfit <- rmeanglm(Surv(time, status) ~ rx, time = 2500, data = colon)
summary(colon.rmfit)
se.rm <- sqrt(diag(vcov(colon.rmfit, type = "robust")))
b.rm <- colon.rmfit$coefficients
conf.rm <- confint(colon.rmfit)

## -----------------------------------------------------------------------------
cbind(eventglm = b.rm, 
      survfit = c(colon.sfit$table[1, 5], 
colon.sfit$table[2:3, 5] - colon.sfit$table[1, 5]))


## -----------------------------------------------------------------------------
colon.ci.adj <- cumincglm(Surv(time, status) ~ rx + age + node4, time = 2500, data = colon)
colon.rm.adj <- rmeanglm(Surv(time, status) ~ rx + age + node4, time = 2500, data = colon)
summary(colon.rm.adj)

## ---- eval = 2----------------------------------------------------------------
?mgus2
head(mgus2)

## -----------------------------------------------------------------------------
crfit <- survfit(Surv(etime, event) ~ sex, eventglm::mgus2)
summary(crfit, times = 120)
print(crfit, rmean = 120)

plot(crfit, col=1:2,  noplot="",
     lty=c(3,3,2,2,1,1), lwd=2, xscale=12,
     xlab="Years post diagnosis", ylab="P(state)")
legend(240, .65, c("Female, death", "Male, death", "malignancy", "(s0)"),
       lty=c(1,1,2,3), col=c(1,2,1,1), bty='n', lwd=2)
abline(v = 120, lty = 2)

## -----------------------------------------------------------------------------
mgfitci <- cumincglm(Surv(etime, event) ~ sex, cause = "pcm", time = 120, 
                   data = mgus2)
summary(mgfitci)

mgfitrmean <- rmeanglm(Surv(etime, event) ~ sex, cause = "pcm", time = 120, 
                       data = mgus2)
summary(mgfitrmean)

## -----------------------------------------------------------------------------
mgfitci2 <- cumincglm(Surv(etime, event) ~ sex + age + hgb, cause = "pcm", 
                      time = 120, data = mgus2)
mgfitrmean2 <- rmeanglm(Surv(etime, event) ~ sex + age + hgb, cause = "pcm", 
                      time = 120, data = mgus2)
summary(mgfitrmean2)

## -----------------------------------------------------------------------------
nboot <- 100 # use a bigger number for real
bootests <- matrix(NA, nrow = nboot, ncol = 4)
for(i in 1:nboot) {
  mgus.b <- mgus2[sample(1:nrow(mgus2), replace = TRUE), ]
  mgfitrmean.b <- rmeanglm(Surv(etime, event) ~ sex + age + hgb, cause = "pcm", 
                      time = 120, data = mgus.b)
  bootests[i,] <- mgfitrmean.b$coefficients
}

se.boot <- sqrt(diag(cov(bootests)))
knitr::kable(cbind(se.boot = se.boot, 
      se.robust = sqrt(diag(vcov(mgfitrmean2))), 
      #se.corrected = sqrt(diag(vcov(mgfitrmean2, type = "corrected"))), 
      se.naive = sqrt(diag(vcov(mgfitrmean2, type = "naive")))), digits = 3)


## -----------------------------------------------------------------------------
hist(predict(mgfitrmean2, newdata = mgus2), 
     xlab = "Predicted lifetime lost due to PCM", main = "")

mgus2$prob.pcm10 <- predict(mgfitci2, newdata = mgus2)
mgus2$pseudo.ci <- mgfitci$y
summary(mgus2$prob.pcm10)
cutps <- quantile(mgus2$prob.pcm10, seq(.1, .9, by = .1), na.rm = TRUE)
mgus2$prob.cut <- cut(mgus2$prob.pcm10, 
                      cutps)

pred.p <- cutps[-length(cutps)] + diff(cutps)
obs.p <- c(by(mgus2$pseudo.ci, mgus2$prob.cut, mean))

plot(obs.p ~ pred.p, xlab = "predicted", ylab = "observed")
abline(0, 1)

## -----------------------------------------------------------------------------
colon.ci.cen1 <- cumincglm(Surv(time, status) ~ rx + age + node4, time = 2500, 
                           data = colon, model.censoring = "stratified", 
                           formula.censoring = ~ rx)

## -----------------------------------------------------------------------------
colon.ci.cen2 <- cumincglm(Surv(time, status) ~ rx + age + node4, time = 2500, 
                           data = colon, model.censoring = "coxph", 
                           formula.censoring = ~ rx + age + node4)
colon.ci.cen3 <- cumincglm(Surv(time, status) ~ rx + age + node4, time = 2500, 
                           data = colon, model.censoring = "aareg", 
                           formula.censoring = ~ rx + age + node4)

knitr::kable(cbind("indep" = colon.ci.adj$coefficients, 
      "strat" = colon.ci.cen1$coefficients, 
      "coxipcw" = colon.ci.cen2$coefficients, 
      "aalenipcw" = colon.ci.cen3$coefficients), digits = 3)

## -----------------------------------------------------------------------------
library(data.table)

# from https://pclambert.net/data/colon.dta
colon2 <- rio::import(system.file("extdata", "colon.dta", package = "eventglm"))
colon2$surv_mm_trunc <- ifelse(colon2$surv_mm > 120.5, 120.5, colon2$surv_mm)

colon2$death <- colon2$status %in% c(1, 2)
colon2$death[colon2$surv_mm > colon2$surv_mm_trunc] <- 0

# from https://pclambert.net/data/popmort.dta
lifetab <- data.table(rio::import(system.file("extdata", "popmort.dta", package = "eventglm")))

## -----------------------------------------------------------------------------

lifetab[, prob.5 := prod(lifetab[`_age` %in% .BY[["_age"]]:(.BY[["_age"]]+4) & 
                              `_year` %in% .BY[["_year"]]:(.BY[["_year"]]+4) & 
                                `_year` == .BY[["_year"]] - .BY[["_age"]] + `_age` &
                              sex == .BY[["sex"]] ]$prob, na.rm = TRUE), 
        by = c("sex", "_year", "_age")]

lifetab[, prob.10 := prod(lifetab[`_age` %in% .BY[["_age"]]:(.BY[["_age"]]+9) & 
                              `_year` %in% .BY[["_year"]]:(.BY[["_year"]]+9) & 
                                `_year` == .BY[["_year"]] - .BY[["_age"]] + `_age` &
                              sex == .BY[["sex"]] ]$prob, na.rm = TRUE), 
        by = c("sex", "_year", "_age")]


colon2 <- merge(colon2, lifetab, 
                by.x = c("sex", "yydx", "age"), 
                by.y = c("sex",  "_year", "_age"), all.x = TRUE, all.y = FALSE)

fit1 <- cumincglm(survival::Surv(surv_mm_trunc, death) ~ 1, data = colon2, time = 1 * 12)
fit5 <- cumincglm(survival::Surv(surv_mm_trunc, death) ~ 1, data = colon2, time = 5 * 12)
fit10 <- cumincglm(survival::Surv(surv_mm_trunc, death) ~ 1, data = colon2, time = 10 * 12)

colon2$po_1 <- 1 - fit1$y
colon2$po_5 <- 1 - fit5$y
colon2$po_10 <- 1 - fit10$y

knitr::kable(cbind(time = c(1, 5, 10), 
                   relsurv.pseudo = with(colon2, 
                                  c(mean(po_1 / prob), 
                                    mean(po_5 / prob.5), 
                                    mean(po_10 / prob.10)),
                                  ),
                   relsurv.pohar = c(0.682, 0.479, 0.441)), digits = 3)

## -----------------------------------------------------------------------------
colon2$status2 <- factor(ifelse(colon2$status == 4, 0, colon2$status), 
                         labels = c("censored", "cancer death", "other death"))
subc <- rbinom(nrow(colon2), size = 1, p = .2)
samp.ind <- subc + (1 - subc) * (colon2$status == 1) * rbinom(nrow(colon2), size = 1, p = .9)
colon.cc <- colon2[as.logical(samp.ind), ]
colon.cc$samp.wt <- 1 / ifelse(colon.cc$status == 1, .2 + .8 * .9, .2)

## -----------------------------------------------------------------------------
cfit.cc <- cumincglm(Surv(surv_mm, status2) ~ age + sex + factor(subsite), 
                     cause = "cancer death", time = 5 * 12, data = colon.cc, 
                     weights = samp.wt)
cfit.full <- cumincglm(Surv(surv_mm, status2) ~ age + sex + factor(subsite), 
                     cause = "cancer death", time = 5 * 12, data = colon2)
knitr::kable(cbind(casecohort = cfit.cc$coefficients, 
      fullsamp = cfit.full$coefficients), digits = 2)

