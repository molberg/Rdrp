## align

#' Examine a dataset or single spectrum
#'
#' Print class and dimension information.
#' @param s a dataset (class 'spectra') or a single spectrum (class 'spectrum')
#' @export
examine <- function(s) {
    print(names(s$head))
    print(s$head)
    print(class(s$head))
    print(class(s$freq))
    if (class(s$freq) == "matrix") print(dim(s$freq))
    else                           print(length(s$freq))
    ## print(head(s$freq))
    print(class(s$data))
    if (class(s$data) == "matrix") print(dim(s$data))
    else                           print(length(s$data))
    ## print(head(s$data))
}

#' Apply function to each spectrum in a dataaset
#'
#' Given a dataset of class 'spectra' and a function, apply this
#' function to each of the spectra contained in the dataset.
#' @param L a list of spectra
#' @param FUN a function operating on a single spectrum
#' @param ... further arguments to FUN
#' @param simplify if TRUE, try to unlist the result and return as matrix or vector.
#' @export
drapply <- function(L, FUN, ..., simplify=TRUE) {
    if (!class(L) == "spectra") {
        warning("first argument should be of class 'spectra''")
        return (NULL)
    }
    result <- list()
    ns <- length(L)
    for (i in seq(ns)) {
        S <- L[[i]]
        result[[i]] <- FUN(S, ...)
    }
    ## print(names(result))
    ## print(class(result[[1]]))
    if (simplify) {
        result <- unlist(result)
        dims <- getDimension(L)
        ## print(dims)
        if (length(result) == prod(dims)) {
            result <- matrix(result, nrow=dims[1], ncol=dims[2])
        }
    }
    result
}

#' Fit a baseline to spectra in a dataset
#'
#' For all spectra in a dataset fit a polynomial of given order, using only the
#' spectral channels given in window.
#' @param s a single spectrum
#' @param order the order of the polynomial to fit
#' @param mask a channel line mask, fit baseline to channels where mask = FALSE
#' @return the fitted baselines
#' @export
baseline <- function(s, order=1, mask=NULL) {
    x <- seq(-1, 1, length.out=length(s$data))
    if (!is.null(mask)) {
        xb <- x[!mask]
        yb <- s$data[!mask]
    } else {
        xb <- x
        yb <- s$data
    }
    fit <- lm(yb ~ poly(xb, order, raw=TRUE))
    res <- as.numeric(residuals(fit))
    print(sd(res, na.rm=TRUE))
    bl <- as.numeric(predict(fit, newdata=data.frame(xb=x)))
    bl
}

## fold
## redres

#' Plot a set of positions
#'
#' Given two vectors lon and lat, which describe the longitudinal and
#' latitudinal coordinates of individual spectra, produce a scatter plot
#' taking into account projection effects.
#'
#' @param lon a numeric vector with longitudinal values, e.g. RA
#' @param lat a numeric vector with latitudinal values, e.g. Dec
#' @param ... further arguments to be passed on to generic function plot
positions <- function(lon, lat, ...) {
    proj <- cos(mean(lat, na.rm=TRUE)*pi/180)
    plot(lon, lat, type='p', asp=1/proj, ...)
    grid()
}

#' Plot a grid of spectra
#'
#' Try to find a grid of rectangles large enough to hold all the
#' spectra and as square like as possible, and "stamp" the individual
#' spectra into the grid cells.
#'
#' @param L a list of spectra
#' @export
stamp <- function(L) {
    op <- par(mar=c(0.0,0.0,0.0,0.0), oma=c(3,3,3,3))
    ns <- length(L)
    nx <- floor(sqrt(ns))
    ny <- nx
    while (nx*ny < ns) {
        nx <- nx+1
    }
    freq <- getFreq(L)
    data <- getData(L)

    par(mfrow=c(ny, nx))
    ymin=min(data, na.rm=TRUE)
    ymax=max(data, na.rm=TRUE)
    for (i in seq(ns)) {
        plot(freq[,i], data[,i], type='l', col='blue',
             ylim=c(ymin,ymax), xaxt='n', yaxt='n', xlab="", ylab="")
    }
    par(mfrow=c(1, 1))
    par(op)
}
