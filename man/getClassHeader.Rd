% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getClassHeader}
\alias{getClassHeader}
\title{Get header infromation from a GILDAS/CLASS single dish data file}
\usage{
getClassHeader(filename)
}
\arguments{
\item{filename}{name of the GILDAS file including path to be opened}
}
\value{
data frame with n rows, where n is the number of scans found
}
\description{
Given a filename, open the file and scan it for single dish spectra
or continuum scans. A data frame is returned where each row corresponds
to the header information of one scan.
}
