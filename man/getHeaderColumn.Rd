% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getHeaderColumn}
\alias{getHeaderColumn}
\title{Get one column of header data.}
\usage{
getHeaderColumn(L, index)
}
\arguments{
\item{L}{a list of spectra, each with components 'head', 'freq' and 'data'}

\item{index}{the index of the column to be retrieved}
}
\value{
a data.frame formed by row-binding all the individual 'head's
}
\description{
Given a list L, where each list member is itself a list with
components 'head' (which is a list or a dataframe with one row),
'freq' (a numeric vector) and 'data' (another numeric vector of same
length as freq), return all the 'head' components as a data frame,
which will have as many rows as the length of the original list.
}
\details{
It is assumed that all 'head' components have the same number of
members with identical names and types.
}
\examples{
S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)),
           freq=-5:5, data=rnorm(11))
S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)),
           freq=-5:5, data=rnorm(11))
L <- list(S1,S2)
class(L) <- "spectra"
i <- which(names(S1$head) == "ra")
getHeaderColumn(L, i)

# will result in

# [1] 1.23 5.43

}
\seealso{
\code{\link{modify}}
}
