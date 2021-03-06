% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_ipcw_pos}
\alias{calc_ipcw_pos}
\title{Compute inverse probability of censoring weights pseudo observations}
\usage{
calc_ipcw_pos(mr, time, causen, type, ipcw.method, Gi)
}
\arguments{
\item{mr}{Model response object returned by \link{Surv}}

\item{time}{Max time}

\item{causen}{Cause of interest (numeric)}

\item{type}{Outcome type, "cuminc", "survival", or "rmean"}

\item{ipcw.method}{"binder" or "hajek"}

\item{Gi}{vector of estimated censoring probabilities}
}
\description{
Compute inverse probability of censoring weights pseudo observations
}
