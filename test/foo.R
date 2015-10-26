library(Rcpp)
library(rbenchmark)

fake <- function(id) {
    head <- data.frame(id=id, RA=0.0, Dec=0.0, T.sys=300.0+id, dt=10.0+id)
    freq <- seq(-500,500)
    data <- rnorm(length(freq))
    list(head=head, freq=freq, data=data)
}

average <- function(L) {
    w <- sapply(L, function(S) { S$head$dt/S$head$T.sys^2 })
    w <- w/sum(w)
    data <- lapply(L, function(S) { as.numeric(S$data) })
    mat <- do.call(cbind, data)
    A <- apply(scale(mat, center=FALSE, scale=1/w), 1, sum)
    A
}

average2 <- function(L) {
    w <- sapply(L, function(s) { s$head$dt/s$head$T.sys^2 })
    w <- w/sum(w)
    data <- lapply(L, function(s) { as.numeric(s$data) })
    mat <- do.call(cbind, data)
    y <- accum(mat, w)
    y
}


target <- "L183"
run <- "RUN2"
cmd <- sprintf("find /precise/data/O2014a-07/%s -name '*s.fits' -exec grep -l %s {} \\;", run, target)
print(cmd)
fitsfiles = system(cmd, intern=TRUE)
L <- Rdrp::readOSO20m(fitsfiles)

benchmark(x <- average(L), y <- average2(L))
f <- L[[1]]$freq
d <- L[[1]]$data
plot(f, d, type='s', col='blue')
lines(f, x, col='red')
lines(f, y, col='green')
lines(f, y-x, type='s', col='black')
