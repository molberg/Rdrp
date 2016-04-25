#library(Rcpp)
#library(rbenchmark)

gauss <- function() {
    head <- data.frame(id=1234, target='foo', line=NA, RA=0.0, Dec=0.0, T.sys=300.0, dt=10.0, f0 = 1421.0)
    freq <-  -200:200*0.01
    data <-  exp(-(freq^2/0.5))+rnorm(length(freq), sd=0.02)
    freq <- freq+1421.0
    G <- list(head=head, freq=freq, data=data)
    class(G) <- "spectrum"
    G
}

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

L <- list()
for (i in seq(20)) {
    L[[i]] <- fake(i)
}

G <- gauss()
plot(G, type='s') #, xlim=c(1420.5,1421.5))
F <- seq(min(G$freq), max(G$freq), length.out = length(G$freq)-199)
H = resample(G, F)
points(H$freq, H$data, type='b', col='red')
