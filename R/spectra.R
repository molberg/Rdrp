#' Print method for class spectra.
#' 
#' Simply print the header (data.frame) part of the spectra.
#' @param ds a dataset
#' @export
print.spectra <- function(ds, ...) {
    print(ds$head)
}

#' Plot method for class spectra.
#' 
#' Plot all spectra in a dataset using matplot.
#' @param ds a dataset
#' @export
plot.spectra <- function(ds, index=seq(ncol(ds$data)), grid=NULL) {
    if (is.null(grid)) {
        x <- ds$freq[,index]
        y <- ds$data[,index]
        matplot(x, y, type='l', lty=1, xlab="", ylab="")
    } else {
        oldpar <- par(mfrow=grid, mar=c(2,2,1,1))
        for (i in index) {
            plot(ds$freq[,i], ds$data[,i], type='l', col='blue', 
                 xlab="", ylab="")
            mtext(ds$head$target[i], 3, -2, adj=0.02)
        }
        par(oldpar)
    }
}

summary.spectra <- function(object, ...) {
    list(header=summary(object$head), data=summary(object$data))
}

positions <- function(sd, pch=1, col='blue') {
    RA <- sd$head$RA
    Dec <- sd$head$Dec
    proj <- cos(mean(sd$head$Dec, na.rm=TRUE)*pi/180)
    plot(RA, Dec, type='p', pch=pch, col=col, asp=1/proj)
    grid()
}

stamp <- function(sd) {
    op <- par(mar=c(0.0,0.0,0.0,0.0), oma=c(3,3,3,3))
    ns <- nrow(sd$head)
    nx <- floor(sqrt(ns))
    ny <- nx
    while (nx*ny < ns) {
        nx <- nx+1
    }
    par(mfrow=c(ny, nx))
    ymin=min(sd$data, na.rm=TRUE)
    ymax=max(sd$data, na.rm=TRUE)
    for (i in seq(ns)) {
        plot(sd$freq[,i], sd$data[,i], type='l', col='blue',
             ylim=c(ymin,ymax), xaxt='n', yaxt='n', xlab="", ylab="")
    }
    par(mfrow=c(1, 1))
    par(op)
}
