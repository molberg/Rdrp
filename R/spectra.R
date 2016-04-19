#' Print method for class spectra.
#'
#' Simply print the header (data.frame) part of the spectra.
#' @param L a list of spectra
#' @export
print.spectra <- function(L, ...) {
    print(getHead(L))
}

#' Plot method for class spectra.
#'
#' Plot all spectra in a dataset using matplot or a grid of plots.
#' @param L a list of spectra
#' @export
plot.spectra <- function(L, index=seq(ncol(ds$data)), grid=NULL) {
    head <- getHead(L)
    x <- getFreq(L)
    y <- getData(L)
    if (drp.options$system == "velocity") {
        x <- velocity(x, head$f0, head$v.LSR)
    }
    if (is.null(grid)) {
        matplot(x[,index], y[,index], type='l', lty=1, xlab="", ylab="")
    } else {
        oldpar <- par(mfrow=grid, mar=c(2,2,1,1))
        for (i in index) {
            plot(x[,i], y[,i], type='l', col='blue', xlab="", ylab="")
            if ("mask" %in% names(ds)) {
                j <- ds$mask[,i]
                lines(x[j,i], y[j,i], col='green')
            }
            mtext(head$target[i], 3, -2, adj=0.02)
        }
        par(oldpar)
    }
}

positions <- function(lon, lat, pch=1, col='blue') {
    proj <- cos(mean(lat, na.rm=TRUE)*pi/180)
    plot(lon, lat, type='p', pch=pch, col=col, asp=1/proj)
    grid()
}

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
