% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/just_cv.r
\name{just_cv}
\alias{just_cv}
\title{Do cross-validation}
\usage{
just_cv(
  dat,
  dirs,
  variants = "ALL",
  sampling,
  start_iter,
  expect = "both",
  saveCVs = TRUE
)
}
\arguments{
\item{dat}{Data}

\item{dirs}{Directories}

\item{variants}{Selection of model variants}

\item{sampling}{Sampling settings}

\item{saveCVs}{default TRUE}
}
\description{
Loads existing model objects and does cross-validation. If do_cv = TRUE, returns the performance measures.
}
