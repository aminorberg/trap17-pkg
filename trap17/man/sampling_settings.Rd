% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_settings.r
\name{sampling_settings}
\alias{sampling_settings}
\title{Sampling settings}
\usage{
sampling_settings(
  totsamp = 100,
  trans = 50,
  thn = 1,
  nfolds = 2,
  nchains = 2,
  mod_rl_priors = NA
)
}
\arguments{
\item{totsamp}{Number of MCMC samples}

\item{trans}{Number of transition MCMC samples}

\item{thn}{Thinning of MCMC samples}

\item{nfolds}{Number of folds for cross-validation}
}
\value{
List of settings to be passed on for MCMC sampling
}
\description{
Set settings for MCMC sampling with Hmsc
}
