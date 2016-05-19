library(Rdrp)
context("DRP")

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

test_that("sieve is working", {
    S <- list(head=list(target="foo"), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    S <- Rdrp::sieve(S, rep(1,3))
    expect_equal(S$freq, c(2,4,6,8))
    expect_equal(S$data, c(2,4,6,8))
})

test_that("trim is working", {
    S <- list(head=list(target="foo"), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    keep <- 3:8
    S <- Rdrp::trim(S, keep)
    expect_equal(S$freq, seq(3,8))
    expect_equal(S$data, seq(3,8))
})

test_that("reverse is working", {
    S <- list(head=list(target="foo"), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    S <- Rdrp::reverse(S)
    expect_equal(S$freq, seq(10,1))
    expect_equal(S$data, seq(10,1))
})

test_that("average is working", {
    S1 <- list(target="foo", head=list(dt=10, T.sys=300), freq=seq(10),data=rep(1,10))
    S2 <- list(target="foo", head=list(dt=20, T.sys=300), freq=seq(10),data=rep(2,10))
    class(S1) <- "spectrum"
    class(S2) <- "spectrum"
    L <- list(S1, S2)
    A <- average(L)
    expect_equal(A$data, rep(5/3, 10))
    S1 <- list(head=list(dt=10, T.sys=300), freq=seq(10),data=rep(1,10))
    S2 <- list(head=list(dt=10, T.sys=100), freq=seq(10),data=rep(2,10))
    class(S1) <- "spectrum"
    class(S2) <- "spectrum"
    L <- list(S1, S2)
    A <- average(L)
    expect_equal(A$data, rep(1.9, 10))
})

test_that("resample is working", {
    S <- list(head=list(target="foo", dt=10, T.sys=300), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    f1 <- seq(2.1,9.1,by=1)
    S <- Rdrp::resample(S, f1, smooth=FALSE)
    expect_equal(S$freq, f1)
    expect_equal(S$data, f1)
    S <- list(head=list(target="foo", dt=10, T.sys=300), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    f1 <- seq(0.1,7.1,by=1)
    S <- Rdrp::resample(S, f1, smooth=FALSE)
    expect_equal(S$freq, f1)
    expect_equal(S$data, c(NA, f1[2:8]))
    if (FALSE) {
        f <- seq(182500,183500)
        f0 = 183000
        S <- list(head=list(id=1, target="foo", dt=10, T.sys=300, f0=f0),
                  freq=f, data=rnorm(length(f))+5.0*exp(-(f-f0)^2/100.0))
        class(S) <- "spectrum"
        plot(S, type='s', xlim=c(182800,183200))
        N <- 100
        f1 <- seq(182500, 183500, length.out=N)
        abline(v=f1, lty=2, col='grey')
        ## points(f1, rep(0, N), col='red', pch=3)
        S1 <- resample(S, f1)
        lines(S1$freq, S1$data, lwd=1, col='green')
        points(S1$freq, S1$data, col='green')
        S2 <- resample(S, f1, smooth=FALSE)
        lines(S2$freq, S2$data, lwd=1, col='red')
        points(S2$freq, S2$data, col='red')
        r.df <- sqrt((S1$freq[2]-S1$freq[1])/(S$freq[2]-S$freq[1]))
        r.rms <- sd(S$data)/sd(S1$data)
        percent <- (r.rms-r.df)/r.rms
        print(c(r.df, r.rms, percent))
        expect_lt(percent, 0.20)
    }
})
