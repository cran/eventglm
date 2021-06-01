## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(survival)
library(eventglm)

## -----------------------------------------------------------------------------
eventglm::pseudo_independent

## -----------------------------------------------------------------------------
pseudo_parametric <- function(formula, time, cause = 1, data,
                        type = c("cuminc", "survival", "rmean"),
                        formula.censoring = NULL, ipcw.method = NULL){

  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))
  
  marginal.estimate <- survival::survreg(margformula, data = data, 
                                         dist = "weibull")

  theta <- pweibull(time, shape = 1 / marginal.estimate$scale, 
                    scale = exp(marginal.estimate$coefficients[1]))
  
  theta.i <- sapply(1:nrow(data), function(i) {
    
    me <- survival::survreg(margformula, data = data[-i, ], dist = "weibull")
    pweibull(time, shape = 1 / me$scale, 
                    scale = exp(me$coefficients[1]))
  
    
  })
  
  POi <- theta  + (nrow(data) - 1) * (theta - theta.i)
  POi

}


## -----------------------------------------------------------------------------
fitpara <- cumincglm(Surv(time, status) ~ rx + sex + age, time = 2500, 
                  model.censoring = pseudo_parametric, 
                  data = colon)

fitdef <- cumincglm(Surv(time, status) ~ rx + sex + age, time = 2500, 
                  model.censoring = "independent", 
                  data = colon)

knitr::kable(sapply(list(parametric = fitpara, default = fitdef), 
       coefficients))

## ---- eval = FALSE------------------------------------------------------------
#  fit1 <- cumincglm(Surv(time, status) ~ rx + sex + age, time = 2500,
#                    model.censoring = "parametric",
#                    data = colon)

## -----------------------------------------------------------------------------
pseudo_infjack <- function(formula, time, cause = 1, data,
                        type = c("cuminc", "survival", "rmean"),
                        formula.censoring = NULL, ipcw.method = NULL) {
  
  marginal.estimate2 <- survival::survfit(update.formula(formula, . ~ 1),
                                             data = data, influence = TRUE)

     tdex <- sapply(time, function(x) max(which(marginal.estimate2$time <= x)))

     pstate <- marginal.estimate2$surv[tdex]


     ## S(t) + (n)[S(t) -S_{-i}(t)]
     POi <- matrix(pstate, nrow = marginal.estimate2$n, ncol = length(time), byrow = TRUE) +
         (marginal.estimate2$n) *
         (marginal.estimate2$influence.surv[, tdex + 1])
     
     POi

}

## -----------------------------------------------------------------------------
fitinf <- cumincglm(Surv(time, status) ~ rx + sex + age, time = 2500, 
                  model.censoring = "infjack", 
                  data = colon)

fitdefsurv <- cumincglm(Surv(time, status) ~ rx + sex + age, time = 2500, 
                  survival = TRUE,
                  data = colon)


knitr::kable(sapply(list(infjack = fitinf, default = fitdefsurv), 
       coefficients))

