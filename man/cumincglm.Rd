% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event-model-fit.R
\name{cumincglm}
\alias{cumincglm}
\title{Generalized linear models for cumulative incidence}
\usage{
cumincglm(
  formula,
  time,
  cause = 1,
  link = "identity",
  model.censoring = "independent",
  formula.censoring = NULL,
  ipcw.method = "binder",
  data,
  weights,
  subset,
  na.action,
  offset,
  control = list(...),
  model = FALSE,
  x = TRUE,
  y = TRUE,
  singular.ok = TRUE,
  contrasts = NULL,
  ...
)
}
\arguments{
\item{formula}{A formula specifying the model. The left hand side must be a
\link[survival]{Surv} object specifying a right censored survival or
competing risks outcome. The status indicator, normally 0=alive, 1=dead.
Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
risks, the event variable will be a factor, whose first level is treated as
censoring. The right hand side is the usual linear combination of
covariates.}

\item{time}{Numeric constant specifying the time at which the cumulative
incidence or survival probability effect estimates are desired.}

\item{cause}{Numeric or character constant specifying the cause indicator of
interest.}

\item{link}{Link function for the cumulative incidence regression model.}

\item{model.censoring}{Type of model for the censoring distribution. Options
are "stratified", which computes the pseudo-observations stratified on a
set of categorical covariates, "aareg" for Aalen's additive hazards model,
and "coxph" for Cox's proportional hazards model. With those options, we
assume that the time to event and event indicator are conditionally
independent of the censoring time, and that the censoring model is
correctly specified. If "independent", we assume completely independent
censoring, i.e., that the time to event and covariates are independent of
the censoring time. the censoring time is independent of the covariates in
the model.}

\item{formula.censoring}{A one sided formula (e.g., \code{~ x1 + x2})
specifying the model for the censoring distribution. If NULL, uses the same
mean model as for the outcome.}

\item{ipcw.method}{Which method to use for calculation of inverse
probability of censoring weighted pseudo observations. "binder" the
default, uses the number of observations as the denominator, while the
"hajek" method uses the sum of the weights as the denominator.}

\item{data}{Data frame in which all variables of formula can be interpreted.}

\item{weights}{an optional vector of 'prior weights' to be used in the
fitting process. Should be NULL or a numeric vector.}

\item{subset}{an optional vector specifying a subset of observations to be
used in the fitting process.}

\item{na.action}{a function which indicates what should happen when the data
contain \code{NA}s. The default is set by the \code{na.action} setting of
\link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
'factory-fresh' default is \link[stats]{na.omit}. Another possible value is
NULL, no action. Value \link[stats]{na.exclude} can be useful.}

\item{offset}{this can be used to specify an a priori known component to be
included in the linear predictor during fitting. This should be NULL or a
numeric vector of length equal to the number of cases. One or more
\link[stats]{offset} terms can be included in the formula instead or as
well, and if more than one is specified their sum is used. See
\link[stats]{model.offset}.}

\item{control}{a list of parameters for controlling the fitting process. This
is passed to \link[stats]{glm.control}.}

\item{model}{a logical value indicating whether model frame should be
included as a component of the returned value.}

\item{x}{logical value indicating whether the model matrix used in the
fitting process should be returned as components of the returned value.}

\item{y}{logical value indicating whether the response vector
(pseudo-observations) used in the fitting process should be returned as
components of the returned value.}

\item{singular.ok}{logical; if FALSE a singular fit is an error.}

\item{contrasts}{an optional list. See the contrasts.arg of
\link[stats]{model.matrix.default}.}

\item{...}{Other arguments passed to \link[stats]{glm.fit}}
}
\value{
A pseudoglm object, with its own methods for print, summary, and
  vcov. It inherits from glm, so predict and other glm methods are supported.
}
\description{
Using pseudo observations for the cumulative incidence, this function then
runs a generalized linear model and estimates the parameters representing
contrasts in the cumulative incidence at a particular time (specified by the
\code{time} argument) across covariate values. The link function can be
"identity" for estimating differences in the cumulative incidence, "log" for
estimating ratios, and any of the other link functions supported by
\link[stats]{quasi}.
}
\examples{
    cumincipcw <- cumincglm(Surv(etime, event) ~ age + sex,
         time = 200, cause = "pcm", link = "identity",
         model.censoring = "independent", data = mgus2)
# stratified on only the categorical covariate
     cumincipcw2 <- cumincglm(Surv(etime, event) ~ age + sex,
                         time = 200, cause = "pcm", link = "identity",
                         model.censoring = "stratified",
                         formula.censoring = ~ sex, data = mgus2)
}