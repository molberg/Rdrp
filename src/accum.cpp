#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector accum(NumericMatrix D, NumericVector w) {
    int ncols = D.ncol();
    int nrows = D.nrow();
    NumericVector x(nrows);
    for (int i=0; i < ncols; i++) {
        for (int j=0; j < nrows; j++) {
            x[j] += D[j+i*nrows]*w[i];
        }
    }
    return x;
}
