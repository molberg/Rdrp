#' A package for data reduction of radio astronomical single dish data.
#'
#' The Rdrp package provides commands for the reduction of spectral line
#' data originating from single dish radio astronomical observations. In
#' that sense it is similar to programs like CLASS, which is part of the
#' \url{https://www.iram.fr/IRAMFR/GILDAS/} system. But
#' contrary to CLASS (and other alternatives), Rdrp is built as a
#' package for the R programming language, thus building on a huge
#' reservoir of other numerical and statistical tools. Also, the user
#' gets the high quality, powerful built-in graphics from R for free.
#'
#' The principle design is as follows: data are organized in memory in
#' the form of lists (of class \code{spectrum}), which have three componnents:
#' \describe{
#' \item{\strong{head}}{header information for the spectrum, itself a list}
#' \item{\strong{freq}}{a numeric vector giving the frequencies of all spectral channels}
#' \item{\strong{data}}{a numeric vector giving the data content of all spectral channels}
#' }
#'
#' The available routines to read data from a FITS file or a CLASS
#' formatted binary file, typically return lists of objects, each of
#' which will be of the format described above. Helper functions
#' \code{\link{getHead}}, \code{\link{getFreq}} and
#' \code{\link{getData}} are available to return all headers as a data
#' frame (with one row per spectrum) and all the frequency (or data)
#' vectors as a numeric matrix, where spectral channels run along rows
#' and the columns correspond to an individual spectrum. Such a list of
#' lists would be of class \code{spectra}, i.e. plural of \code{spectrum}.
#'
#' Here is an example of what a short session may look like:
#' @examples
#' library(Rdrp)
#' \dontrun{
#' assign("system","velocity", Rdrp::options) # work in velocity space
#' # now read some data which e.g. we may have obtained at APEX
#' L <- readClass("mydata.apex")
#' H <- getHead(L)
#' print(H)            # take a look at the header information
#' i <- which(H$target == "IC348" & H$line == "CO(3-2)")
#' L <- L[i]           # only keep spectra with given target, line
#' S <- L[[1]]         # get the first spectrum
#' plot(S)             # and plot it; this will use plot.spectrum(...)
#' A <- average(L)     # form the average
#' # we expect a spectral line between -20 and +20 km/s
#' linemask <- mask(A, c(-20,20))
#' # fit a second order baseline
#' bl <- baseline(A, order=2, mask=linemask)
#' plot(A)                           # plot the spectrum ...
#' lines(velocity(A), bl, col='red') # ... and fitted baseline
#' A$data <- A$data - bl             # subtract baseline
#' plot(A, type='s', xlab="velocity [km/s]", ylab=expression(T[A]))
#' }
#' # Here is how you would construct a very simple, fake spectrum:
#' head <- list(target="my target", ra=1.0, dec=2.0, f0=1421.0)
#' freq <- head$f0 + seq(-100,100)*0.1  # +-10 MHz around centre frequency
#' data <- rnorm(length(freq))
#' S <- list(head=head, freq=freq, data=data)
#' class(S) <- "spectrum"
#' plot(S, type='s', col='red')
#'
#' @docType package
#' @name Rdrp-package
#' @author Michael Olberg, \email{michael.olberg@@chalmers.se}
NULL

## NULL
