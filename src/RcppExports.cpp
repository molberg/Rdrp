// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// accum
NumericVector accum(NumericMatrix D, NumericVector w);
RcppExport SEXP Rdrp_accum(SEXP DSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(accum(D, w));
    return __result;
END_RCPP
}
// getScanNumber
int getScanNumber(List ds);
RcppExport SEXP Rdrp_getScanNumber(SEXP dsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type ds(dsSEXP);
    __result = Rcpp::wrap(getScanNumber(ds));
    return __result;
END_RCPP
}
// getZeros
arma::vec getZeros(int n);
RcppExport SEXP Rdrp_getZeros(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(getZeros(n));
    return __result;
END_RCPP
}
