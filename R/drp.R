#getHead <- function(L) {
#    headers <- do.call(rbind, lapply(L, function(x) { as.data.frame(x$head) }))
#    H <- as.data.frame(headers)
#    H
#}

#' Modify a header column
#'
#' Supply new values for a column of the header
#' @param L a list of spectra
#' @param colname a string specifying which header column to modify
#' @param value a vector holding the new values of the header column
#' @export
modify <- function(L, colname, value) {
    if (length(value) == 1) value <- rep(value, length(L))
    if (length(value) != length(L)) stop("wrong length of supplied value")
    L <- lapply(seq(length(L)), function(i) {
        S <- L[[i]]
        S$head[colname] <- value[i]
        S
    })
    L
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
#' @param S a single spectrum
#' @param order the order of the polynomial to fit
#' @param mask a channel line mask, fit baseline to channels where mask = FALSE
#' @return a list giving the fitted baseline and the residuals of the fit
#' @seealso \code{\link{mask}}
#' @examples
#' data(salsa)
#' A <- average(salsa)
#' limits <- c(1419.9, 1421.1)            # define line region
#' linemask <- mask(A, limits)            # turn into mask
#' bl <- baseline(A, order=3, mask=linemask)  # fit baseline
#' plot(A)                                # plot average spectrum
#' abline(v=limits, lty=2, col='grey')    # mark line region
#' print(bl$residuals)                    # print the residuals of the fit
#' lines(A$freq, bl$fitted, col='red')    # draw baseline on top
#' @export
baseline <- function(S, order=1, mask=NULL) {
    x <- seq(-1, 1, length.out=length(S$data))
    if (!is.null(mask)) {
        xb <- x[!mask]
        yb <- S$data[!mask]
    } else {
        xb <- x
        yb <- S$data
    }
    fit <- stats::lm(yb ~ poly(xb, order, raw=TRUE))
    res <- as.numeric(stats::residuals(fit))
    bl <- list(fitted=as.numeric(stats::predict(fit, newdata=data.frame(xb=x))),
               residuals=stats::sd(res, na.rm=TRUE))
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
#' @export
positions <- function(lon, lat, ...) {
    proj <- cos(mean(lat, na.rm=TRUE)*pi/180)
    graphics::plot(lon, lat, type='p', asp=1/proj, ...)
    graphics::grid()
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
    op <- graphics::par(mar=c(0.0,0.0,0.0,0.0), oma=c(3,3,3,3))
    ns <- length(L)
    nx <- floor(sqrt(ns))
    ny <- nx
    while (nx*ny < ns) {
        nx <- nx+1
    }
    freq <- getFreq(L)
    data <- getData(L)

    graphics::par(mfrow=c(ny, nx))
    ymin=min(data, na.rm=TRUE)
    ymax=max(data, na.rm=TRUE)
    for (i in seq(ns)) {
        graphics::plot(freq[,i], data[,i], type='l', col='blue',
                       ylim=c(ymin,ymax), xaxt='n', yaxt='n', xlab="", ylab="")
    }
    graphics::par(mfrow=c(1, 1))
    graphics::par(op)
}

#' Mark a region
#'
#' Use the mouse to mark region(s) of a spectrum which should be avoided
#' in baseline fitting. The number of points to be clicked should be
#' even. The routine will call \code{mask(S, limits)} on exit, with limits
#' being the x-coordinates of all the clicked points.
#'
#' @param S a single spectrum, the one that is currently plotted
#' @param ... further parameters passed to 'identify'
#' @return a logical vector, one per channel
#' @seealso \code{\link{mask}}
#' @export
markRegion <- function(S, ...) {
    x <- S$freq
    if (getOption("system") == "velocity") x <- velocities(S)
    xy <- grDevices::xy.coords(x, S$data)
    x <- xy$x
    n <- length(x)
    y <- xy$y
    sel <- rep(FALSE, n); res <- integer(0)
    while(sum(sel) < n) {
        ans <- graphics::identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
        if (!length(ans)) break
        ans <- which(!sel)[ans]
        ## points(x[ans], y[ans], pch = pch)
        graphics::abline(v=x[ans], lty=2, col='grey')
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
    nw <- length(res)
    if ((nw %% 2) == 1) stop("length of limits must be even")
    W <- matrix(x[res], ncol=2, byrow=TRUE)
    DF <- data.frame(from=W[,1], to=W[,2])
    print(DF)
    usr <- graphics::par("usr")
    y0 <- usr[3]+0.1*(usr[4]-usr[3])
    mx <- c(x[1], rep(x[res], each=2), x[length(x)])
    my <- c(rep(rep(c(0,1), nw/2), each=2), 0, 0)
    graphics::lines(mx, y0+0.1*my*(usr[4]-usr[3]), type='s', col='red')
    M <- mask(S, x[res])
    M
}
