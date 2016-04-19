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
#' @return the fitted baselines
#' @export
baseline <- function(s, order=1) {
    x <- seq(-1, 1, length.out=length(s$data))
    if ("mask" %in% names(s) & !is.null(s[["mask"]])) {
        xb <- x[!s$mask]
        yb <- s$data[!s$mask]
    } else {
        xb <- x
        yb <- s$data
    }
    fit <- lm(yb ~ poly(xb, order, raw=TRUE))
    # print(fit)
    bl <- list(fit=predict(fit, newdata=data.frame(xb=x)), rms=sd(yb, na.rm=TRUE))
    bl
}

#' Apply boxcar smoothing.
#'
#' Filter the data vector by a a boxcar filter function.
#' @param s a single spectrum
#' @param n the length of the boxcar filter
#' @return the filtered data vector
#' @export
boxcar <- function(s, n=3) {
    box <- rep(1,n)/n
    d <- as.double(filter(s$data, box))
    d
}

#' Clip channels from spectra
#'
#' Given a vector of channel numbers, only keep those channels from both the frequency
#' and data matrices.
#' @param ds a datase
#' @param keep an integer vector listing channel numbers to preserve
#' @return a new dataset with fewer channels
#' @export
clipChannels <- function(ds, keep) {
    f <- ds$freq[keep,]
    d <- ds$data[keep,]
    sd <- list(head=ds$head, freq=as.matrix(f), data=as.matrix(d))
    class(sd) <- "spectra"
    sd
}

#' Construct a line mask
#'
#' For a given set of frequency (or velocity) windows, construct a
#' matrix of logical values, defining the areas which should not be used
#' in baseline fitting.
#' @param y a matrix of frequency or velocity vectors, one per column
#' @param limits pairs of values which each define a window
#' @export
lineMask <- function(y, limits) {
    if ((length(limits) %% 2) == 1) {
        warning("length of limits must be even")
        return (NULL)
    }
    window <- matrix(FALSE, nrow=nrow(y), ncol=ncol(y))
    for (i in seq(1, length(limits), by=2)) {
        for (j in seq(ncol(y))) {
            window[,j] <- (y[,j] >= limits[i] & y[,j] <= limits[i+1])
        }
    }
    window
}


## filter
## integrate
## fold
## redres
## smooth

#' Calculate velocity vectors
#'
#' Given the frequency vectors and source velocities, calculate a matrix
#' of velocities corresponding to the ds$freq matrix.
#' @param f a matrix of frequency vectors for which velocities will be calculated
#' @param f0 the rest frequency in the same units as f
#' @param vs the source velocity in km/s
#' @return vel a matrix with one velocity vector per column
#' @export
velocity <- function(f, f0, vs) {
    v <- scale(f, center=f0, scale=-f0/299792.46)
    v <- scale(v, center=-vs, scale=FALSE)
    v
}
