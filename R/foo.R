library(Rcpp)
library(rbenchmark)

cppFunction('
    NumericVector rcppaverage(List L) {
        int n = L.size();
        List l = as<List>(L[0]);
        int m = as<NumericVector>(l["data"]).size();
        NumericVector x(m);
        for (int i=0; i < n; i++) {
            l = as<List>(L[i]);
            x = x + as<NumericVector>(l["data"]);
        }
        for (int j=0; j < m; j++) x[j] /= (double)n;
        return x;
    }')

fake <- function(id) {
    head <- data.frame(id=id, RA=0.0, Dec=0.0)
    freq <- seq(-500,500)
    data <- rnorm(length(freq))
    list(head=head, freq=freq, data=data)
}

average <- function(L) {
    data <- lapply(L,
                   function(S) {
                       as.numeric(S$data)
                   })
    mat <- do.call(cbind, data)
    A <- apply(mat, 1, mean)
    A
}

L <- lapply(seq(100), fake)

#print(x)
benchmark(A <- average(L), x <- rcppaverage(L))
f <- L[[1]]$freq
d <- L[[1]]$data
plot(f, d, type='s', col='blue')
lines(f, A, col='red')
lines(f, x, col='green')
lines(f, A-x, type='s', col='black')
