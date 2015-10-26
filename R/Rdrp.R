## align

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
    fit <- lapply(bl,
                   function(x) {
                       as.numeric(x$fit)
                   })
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
