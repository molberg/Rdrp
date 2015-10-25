#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int getScanNumber(List ds) {
    List head = as<List>(ds["head"]);
    int id = as<int>(head["id"]);
    return id;
}
