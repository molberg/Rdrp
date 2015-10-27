## align

#' Average spectra in a dataset
#'
#' Calculate a weighted sum of all spectra in a dataset.
#' @param ds a dataset
#' @return a new dataset with the averaged spectrum as the only entry
#' @export
average <- function(ds) {
    # perform checks
    tol <- drp.options$position.tolerance # 1.0 arcsec
    if (sd(ds$head$RA) > tol | sd(ds$head$Dec) > tol) {
        warning("position mismatch")
    }
    tol <- drp.options$frequency.tolerance/1.0e6 # MHz
    if (sd(apply(ds$freq, 1, sd)) > tol) {
        warning("frequency mismatch")
    }
    w <- ds$head$dt/ds$head$T.sys^2
    w <- w/sum(w)
    # new header has integration time normalized to 300 K Tsys
    h <- ds$head[1,]
    h$dt <- sum(ds$head$dt*(300/ds$head$T.sys)^2)
    h$T.sys <- 300.0
    names(h)[names(h) == "observed.date"] <- "reduced.date"
    h$reduced.date <- Sys.time()
    # use the first frequency vector
    f <- ds$freq[,1]
    # now calculate average
    # d <- apply(scale(ds$data, center=FALSE, scale=1/w), 1, sum)
    # Rcpp version
    d <- accum(ds$data, w)
    sd <- list(head=h, freq=as.matrix(f), data=as.matrix(d))
    class(sd) <- "spectra"
    sd
}

#' Fit a baseline to spectra in a dataset
#'
#' For all spectra in a dataset fit a polynomial of given order, using only the
#' spectral channels given in window.
#' @param ds a dataset
#' @param order the order of the polynomial to fit
#' @param window if given, a logical vector with TRUE for channels to be used
#' @return a list of fitted baselines
#' @export
baseline <- function(ds, order=1, window=NULL) {
    bl <- lapply(seq(nrow(ds$head)), function(i, ds, order, window) {
        x <- seq(-1, 1, length.out=nrow(ds$data))
        if (is.null(window)) {
            print("no window")
            xb <- x
            yb <- ds$data[,i]
        } else {
            xb <- x[window]
            yb <- ds$data[window,i]
        }
        fit <- lm(yb ~ poly(xb, order, raw=TRUE))
        bl <- list(fit=predict(fit, newdata=data.frame(xb=x)),
                   coefficients=fit$coefficients,
                   rms=sd(yb, na.rm=TRUE))
    }, ds, order, window)
    bl
    fit <- lapply(bl, function(x) { as.numeric(x$fit) })
    mat <- do.call(cbind, fit)
    coeffs <- do.call(rbind, lapply(bl, function(x) { x$coefficients }))
    colnames(coeffs) <- paste("c", seq(ncol(coeffs))-1, sep="")
    rms <- sapply(bl, function(x) { x$rms })
    list(fit=as.matrix(mat), coefficients=coeffs, rms=rms)
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

## filter
## integrate
## fold
## redres
## smooth

#' Calculate velocity vectors
#'
#' Given the frequency vectors and source velocities, calculate a matrix
#' of velocities corresponding to the ds$freq matrix.
#' @param ds dataset for which velocities will be calculated
#' @return vel a matrix with one velocity vector per column
#' @export
velocity <- function(ds) {
    f <- ds$freq
    v <- scale(f, center=ds$head$f0, scale=-ds$head$f0/299792.46)
    v <- scale(v, center=-ds$head$v.LSR, scale=FALSE)
    v
}
