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

test_that("filter is working", {
    S <- list(head=list(), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    S <- Rdrp::filter(S, rep(1,3))
    expect_equal(S$freq, c(2,4,6,8))
    expect_equal(S$data, c(2,4,6,8))
})

test_that("trim is working", {
    S <- list(head=list(), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    keep <- 3:8
    S <- Rdrp::trim(S, keep)
    expect_equal(S$freq, seq(3,8))
    expect_equal(S$data, seq(3,8))
})

test_that("reverse is working", {
    S <- list(head=list(), freq=seq(10),data=seq(10))
    class(S) <- "spectrum"
    S <- Rdrp::reverse(S)
    expect_equal(S$freq, seq(10,1))
    expect_equal(S$data, seq(10,1))
})

test_that("average is working", {
    S1 <- list(head=list(dt=10, T.sys=300), freq=seq(10),data=rep(1,10))
    S2 <- list(head=list(dt=20, T.sys=300), freq=seq(10),data=rep(2,10))
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
