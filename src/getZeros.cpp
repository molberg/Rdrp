#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec getZeros(int n)
{
    arma::vec x;
    x.zeros(n);

    return x;
}
