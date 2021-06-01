## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(eventglm)

colon.cifit <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon)
summary(colon.cifit)
confint(colon.cifit)

