#' Print method for class spectrum
#'
#' Simply print the header part of the spectrum.
#' @param x a single spectrum
#' @param ... further arguments to be passed to generic function print
#' @export
print.spectrum <- function(x, ...) {
    print(as.data.frame(x$head))
    print(cbind(x$freq, x$data))
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
#' @param type the type of plot to be used, default is 'l', i.e. line
#' @param col the colour of the plot, default 'blue'
#' @param ... further arguments to be passed to generic plot function
#' @examples
#' \dontrun{
#' plot(S, type='s')   # use 'histo' mode to plot spectrum
#' }
#' @export
plot.spectrum <- function(x, type='l', col='blue', ...) {
    X <- x$freq
    D <- x$data
    if (exists("system", envir=options)) {
        print(options$system)
        if (options$system == "velocity") {
            X = velocity(x)
        }
    }
    # print(cbind(X,D))
    plot(X, D, type=type, col=col, xlab="", ylab="", ...)
    mtext(paste(x$head$id, x$head$target), 3, 1, adj=0.01)
}

#' Plot method for class spectra
#'
#' Plot all spectra in a dataset using matplot or a grid of plots.
#' @param x a list of spectra
#' @param type the type of plot to be used, default is 'l', i.e. line
#' @param col the colour of the plot, default 'blue'
#' @param ... further arguments to be passed to generic plot function
#' @param grid specify grid (nrow x ncol) on which to organize spectra
#' @examples
#' data(salsa)
#' plot(salsa)
#' plot(salsa, grid=c(29,1))
#' @export
plot.spectra <- function(x, type='l', col='blue', ..., grid=NULL) {
    head <- getHead(x)
    freq <- getFreq(x)
    data <- getData(x)
    if (exists("system", envir=options)) {
        print(options$system)
        if (options$system == "velocity") {
            freq = getVelo(x)
        }
    }
    if (is.null(grid)) {
        matplot(freq, data, type=type, lty=1, xlab="", ylab="", ...)
    } else {
        oldpar <- par(mfrow=grid, mar=c(0,0,0,0), oma=c(2,2,2,2))
        for (i in seq(nrow(head))) {
            plot(freq[,i], data[,i], type=type, col=col, xaxt='n', yaxt='n', xlab="", ylab="", ...)
            # mtext(head$target[i], 3, -2, adj=0.02)
        }
        par(oldpar)
    }
}
