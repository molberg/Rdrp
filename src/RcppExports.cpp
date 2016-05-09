// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getHead
DataFrame getHead(List L);
RcppExport SEXP Rdrp_getHead(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(getHead(L));
    return __result;
END_RCPP
}
// getFreq
NumericMatrix getFreq(List L);
RcppExport SEXP Rdrp_getFreq(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(getFreq(L));
    return __result;
END_RCPP
}
// velocity
NumericVector velocity(List S);
RcppExport SEXP Rdrp_velocity(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    __result = Rcpp::wrap(velocity(S));
    return __result;
END_RCPP
}
// getVelo
NumericMatrix getVelo(List L);
RcppExport SEXP Rdrp_getVelo(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(getVelo(L));
    return __result;
END_RCPP
}
// getData
NumericMatrix getData(List L);
RcppExport SEXP Rdrp_getData(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(getData(L));
    return __result;
END_RCPP
}
// getDimension
IntegerVector getDimension(List L);
RcppExport SEXP Rdrp_getDimension(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(getDimension(L));
    return __result;
END_RCPP
}
// average
List average(List L);
RcppExport SEXP Rdrp_average(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    __result = Rcpp::wrap(average(L));
    return __result;
END_RCPP
}
// modify
void modify(List L, std::string column, SEXP value);
RcppExport SEXP Rdrp_modify(SEXP LSEXP, SEXP columnSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< std::string >::type column(columnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type value(valueSEXP);
    modify(L, column, value);
    return R_NilValue;
END_RCPP
}
// addColumns
void addColumns(List L, CharacterVector newnames);
RcppExport SEXP Rdrp_addColumns(SEXP LSEXP, SEXP newnamesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type newnames(newnamesSEXP);
    addColumns(L, newnames);
    return R_NilValue;
END_RCPP
}
// bar
SEXP bar(StringVector which);
RcppExport SEXP Rdrp_bar(SEXP whichSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< StringVector >::type which(whichSEXP);
    __result = Rcpp::wrap(bar(which));
    return __result;
END_RCPP
}
// foo
List foo(List S);
RcppExport SEXP Rdrp_foo(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    __result = Rcpp::wrap(foo(S));
    return __result;
END_RCPP
}
// fold
List fold(List S, double ft, bool shift);
RcppExport SEXP Rdrp_fold(SEXP SSEXP, SEXP ftSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< bool >::type shift(shiftSEXP);
    __result = Rcpp::wrap(fold(S, ft, shift));
    return __result;
END_RCPP
}
// reverse
List reverse(List S);
RcppExport SEXP Rdrp_reverse(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    __result = Rcpp::wrap(reverse(S));
    return __result;
END_RCPP
}
// area
double area(List S, LogicalVector mask);
RcppExport SEXP Rdrp_area(SEXP SSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type mask(maskSEXP);
    __result = Rcpp::wrap(area(S, mask));
    return __result;
END_RCPP
}
// moment
NumericVector moment(List S, LogicalVector mask);
RcppExport SEXP Rdrp_moment(SEXP SSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type mask(maskSEXP);
    __result = Rcpp::wrap(moment(S, mask));
    return __result;
END_RCPP
}
// trim
List trim(List S, IntegerVector keep);
RcppExport SEXP Rdrp_trim(SEXP SSEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type keep(keepSEXP);
    __result = Rcpp::wrap(trim(S, keep));
    return __result;
END_RCPP
}
// filter
List filter(List S, NumericVector coeffs);
RcppExport SEXP Rdrp_filter(SEXP SSEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coeffs(coeffsSEXP);
    __result = Rcpp::wrap(filter(S, coeffs));
    return __result;
END_RCPP
}
// resample
List resample(List S, NumericVector f, bool smooth);
RcppExport SEXP Rdrp_resample(SEXP SSEXP, SEXP fSEXP, SEXP smoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type smooth(smoothSEXP);
    __result = Rcpp::wrap(resample(S, f, smooth));
    return __result;
END_RCPP
}
// rescale
List rescale(List S, double factor, double bias);
RcppExport SEXP Rdrp_rescale(SEXP SSEXP, SEXP factorSEXP, SEXP biasSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< double >::type bias(biasSEXP);
    __result = Rcpp::wrap(rescale(S, factor, bias));
    return __result;
END_RCPP
}
// mask
LogicalVector mask(List S, NumericVector limits);
RcppExport SEXP Rdrp_mask(SEXP SSEXP, SEXP limitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type limits(limitsSEXP);
    __result = Rcpp::wrap(mask(S, limits));
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
// fooData
arma::mat fooData(Rcpp::List L);
RcppExport SEXP Rdrp_fooData(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    __result = Rcpp::wrap(fooData(L));
    return __result;
END_RCPP
}
