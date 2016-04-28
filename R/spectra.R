drp.options <- NULL

#' Print method for class spectrum
#'
#' Simply print the header part of the spectrum.
#' @param x a single spectrum
#' @param ... further arguments to be passed to generic function print
#' @export
print.spectrum <- function(x, ...) {
    print(as.data.frame(x$head))
    print(head(cbind(x$freq, x$data)))
    print(tail(cbind(x$freq, x$data)))
}

#' Print method for class spectra
#'
#' Simply print the header (data.frame) part of the spectra.
#' @param x a list of spectra
#' @param ... further arguments to be passed to generic function print
#' @export
print.spectra <- function(x, ...) {
    print(getHead(x))
}

#' Plot method for class spectrum
#'
#' Plot a single spectrum.
#' @param x a single spectrum
#' @param ... further arguments to be passed to generic function plot
#' @export
plot.spectrum <- function(x, ...) {
    X <- x$freq
    D <- x$data
    if (exists("system", envir=drp.options)) {
        print(drp.options$system)
        if (drp.options$system == "velocity") {
            X = velocity(x)
        }
    }
    plot(X, D, xlab="", ylab="", ...)
    mtext(paste(x$head$id, x$head$target), 3, 1, adj=0.01)
}

#' Plot method for class spectra
#'
#' Plot all spectra in a dataset using matplot or a grid of plots.
#' @param x a list of spectra
#' @param ... further arguments to be passed to generic function plot
#' @param grid specify grid (nrow x ncol) on which to organize spectra
#' @export
plot.spectra <- function(x, ..., grid=NULL) {
    head <- getHead(x)
    freq <- getFreq(x)
    data <- getData(x)
    if (drp.options$system == "velocity") {
        freq = getVelo(x)
    }
    if (is.null(grid)) {
        matplot(freq, data, type='l', lty=1, xlab="", ylab="", ...)
    } else {
        oldpar <- par(mfrow=grid, mar=c(0,0,0,0), oma=c(2,2,2,2))
        for (i in seq(nrow(head))) {
            plot(freq[,i], data[,i], type='l', col='blue', xaxt='n', yaxt='n', xlab="", ylab="", ...)
            # mtext(head$target[i], 3, -2, adj=0.02)
        }
        par(oldpar)
    }
}
