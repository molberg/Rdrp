% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drp.R
\name{drapply}
\alias{drapply}
\title{Apply function to each spectrum in a dataaset}
\usage{
drapply(L, FUN, ..., simplify = TRUE)
}
\arguments{
\item{L}{a list of spectra}

\item{FUN}{a function operating on a single spectrum}

\item{...}{further arguments to FUN}

\item{simplify}{if TRUE, try to unlist the result and return as matrix or vector.}
}
\description{
Given a dataset of class 'spectra' and a function, apply this
function to each of the spectra contained in the dataset.
}
