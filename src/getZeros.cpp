#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec getZeros(int n)
{
    arma::vec x;
    x.zeros(n);

    return x;
}

// [[Rcpp::export]]
arma::mat fooData(Rcpp::List L) 
{
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    Rcpp::List l = L[0];

    /* get the data part, their length is the number of channels */
    arma::vec data = Rcpp::as<arma::vec>(l["data"]);
    int nChannels = data.size();

    /* allocate matrix */
    arma::mat D(nChannels, nSpectra);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        data = Rcpp::as<arma::vec>(l["data"]);
        if (data.size() != nChannels) {
            sprintf(error, "data vector %d has wrong length: %d <> %d",
                    row, data.size(), nChannels);
            Rcpp::stop(error);
        }
        D.col(row) = data;
        // for (int k = 0; k < nChannels; k++) {
        //     D(k,row) = data[k];
        // }
    }

    return D;
}
