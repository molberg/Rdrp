library(Rcpp)
library(rbenchmark)

fake <- function(id) {
    head <- data.frame(id=id, target='foo', line=NA, RA=0.0, Dec=0.0, T.sys=300.0+id, dt=10.0+id)
    freq <- seq(-128,128)
    data <- rnorm(length(freq))
    list(head=head, freq=freq, data=data)
}

average1 <- function(ds) {
    w <- ds$head$dt/ds$head$T.sys^2
    w <- w/sum(w)
    A <- apply(scale(ds$data, center=FALSE, scale=1/w), 1, sum)
    A
}

average2 <- function(ds) {
    w <- ds$head$dt/ds$head$T.sys^2
    w <- w/sum(w)
    y <- accum(ds$data, w)
    y
}


#target <- "L183"
#run <- "RUN2"
#cmd <- sprintf("find /precise/data/O2014a-07/%s -name '*s.fits' -exec grep -l %s {} \\;", run, target)
#print(cmd)
#fitsfiles = system(cmd, intern=TRUE)
#L <- Rdrp::readOSO20m(fitsfiles)
#ds <- createDataset(L)

#benchmark(x <- average1(ds), y <- average2(ds))
#f <- L[[1]]$freq
#d <- L[[1]]$data
#plot(f, d, type='s', col='blue')
#lines(f, x, col='red')
#lines(f, y, col='green')
#lines(f, y-x, type='s', col='black')

#A <- average(ds)
#plot(A)
#print(A)

S <- list()
for (i in seq(20)) {
    S[[i]] <- fake(i)
    plot(S[[i]]$freq, S[[i]]$data, type='s', col='blue')
}
ds <- createDataset(S)
print(ds$head)
print(drapply(ds))
