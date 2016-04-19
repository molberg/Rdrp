#include <Rcpp.h>

using namespace Rcpp;

//' Get header data.
//'
//' Given a list L, where each list member is itself a list with components 'head' (which is a list
//' or a dataframe with one row), 'freq' (a numeric vector) and 'data' (another numeric vector of same
//' length as freq), return a dataframe for the header information, with as many rows as the length
//' of the original list.
//'
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @return a data.frame having all the individual headers as rows
//' @examples
//' S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)), freq=-5:5, data=rnorm(11))
//' S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)), freq=-5:5, data=rnorm(11))
//'
//' getHead(list(S1,S2))
//'
//' will result in
//'
//'   target   ra  dec dt
//' 1  Orion 1.23 -0.5 20
//' 2  SgrB2 5.43  0.5 20
//'
// [[Rcpp::export]]
DataFrame getHead(List L) {
    static char error[80];
    NumericVector d;
    IntegerVector i;
    CharacterVector c;

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the header part, its length is the number of header columns */
    List h = l[0];
    int nHeader = h.length();

    /* get header column names */
    StringVector names = h.names();
    /* vector to hold column types */
    IntegerVector types(nHeader);

    /* create a list which will be casted into a dataframe later */
    List H(nHeader);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        /* get the header part */
        h = l[0];
        /* if this is not the first spectrum we are parsing, make sure we have correct number of columns */
        if (row > 0 && h.length() != nHeader) {
            sprintf(error, "header %d has wrong number of columns: %ld <> %d", row, h.length(), nHeader);
            stop(error);
        }
        for (int col = 0; col < nHeader; col++) {
            SEXP a = h[col];
            int typ = TYPEOF(a);
            switch (typ) {
              case LGLSXP: /* logical header data */
                if (row == 0) {
                    H[col] = LogicalVector(nSpectra);
                    types[col] = LGLSXP;
                } else {
                    if (types[col] != typ) {
                        sprintf(error, "type mismatch in column %d: %d <> %d", col, types[col], typ);
                        stop(error);
                    }
                }
                c = H[col];
                c[row] = as<bool>(h[col]);
                break;
              case INTSXP: /* integer header data */
                /* if this is first spectrum, retrieve type and create empty vector for storing the column */
                if (row == 0) {
                    H[col] = IntegerVector(nSpectra);
                    types[col] = INTSXP;
                } else {
                    if (types[col] != typ) {
                        sprintf(error, "type mismatch in column %d: %d <> %d", col, types[col], typ);
                        stop(error);
                    }
                }
                i = H[col];
                i[row] = as<int>(h[col]);
                break;
              case REALSXP: /* floating point (double) header data */
                if (row == 0) {
                    H[col] = NumericVector(nSpectra);
                    types[col] = REALSXP;
                } else {
                    if (types[col] != typ) {
                        sprintf(error, "type mismatch in column %d: %d <> %d", col, types[col], typ);
                        stop(error);
                    }
                }
                d = H[col];
                d[row] = as<double>(h[col]);
                break;
              case STRSXP: /* string header data */
                if (row == 0) {
                    H[col] = CharacterVector(nSpectra);
                    types[col] = STRSXP;
                } else {
                    if (types[col] != typ) {
                        sprintf(error, "type mismatch in column %d: %d <> %d", col, types[col], typ);
                        stop(error);
                    }
                }
                c = H[col];
                c[row] = as<std::string>(h[col]);
                break;
              case VECSXP: /* generic vector header data */
                if (row == 0) {
                    H[col] = IntegerVector(nSpectra);
                    types[col] = VECSXP;
                } else {
                    if (types[col] != typ) {
                        sprintf(error, "type mismatch in column %d: %d <> %d", col, types[col], typ);
                        stop(error);
                    }
                }
                break;
              default:
                sprintf(error, "invalid column data type: %d", typ);
                warning(error);
                if (row == 0) {
                    H[col] = IntegerVector(nSpectra);
                    types[col] = INTSXP;
                }
                i = H[col];
                i[row] = NA_INTEGER;
                // sprintf(error, "invalid column data type: %d", typ);
                // stop(error);
                break;
            }
        }
    }

    /* So now we have all the headers in list H */
    H.attr("names") = names;
    /* turn H into a dataframe */
    Language call("as.data.frame", H);
    DataFrame DF = call.eval();
    return DF;
}

//' Get frequency vectors
//'
//' From a list of spectra, get the frequency vectors and return as matrix.
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @return a matrix having all the frequency vectors as columns
// [[Rcpp::export]]
NumericMatrix getFreq(List L) {
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector freq = l[1];
    int nChannels = freq.length();

    /* allocate matrix */
    NumericMatrix F(nChannels, nSpectra);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        freq = l[1];
        if (freq.length() != nChannels) {
            sprintf(error, "freq vector %d has wrong length: %ld <> %d",
                    row, freq.length(), nChannels);
            stop(error);
        }
        for (int k = 0; k < nChannels; k++) {
            F(k,row) = freq[k];
        }
    }

    return F;
}

//' Get data vectors
//'
//' From a list of spectra, get the data vectors and return as matrix.
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @return a matrix having all the data vectors as columns
// [[Rcpp::export]]
NumericMatrix getData(List L) {
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector data = l[2];
    int nChannels = data.length();

    /* allocate matrix */
    NumericMatrix D(nChannels, nSpectra);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        data = l[2];
        if (data.length() != nChannels) {
            sprintf(error, "data vector %d has wrong length: %ld <> %d",
                    row, data.length(), nChannels);
            stop(error);
        }
        for (int k = 0; k < nChannels; k++) {
            D(k,row) = data[k];
        }
    }

    return D;
}

//' Get dimensions
//'
//' From a list of spectra, get the dimensions of the data.
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @return a two component integer vector (nChannels, nSpectra)
// [[Rcpp::export]]
IntegerVector getDimension(List L) {
    IntegerVector dim(2);

    /* the length of our list is the number of spectra */
    dim[1] = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector data = l[2];
    dim[0] = data.length();

    return dim;
}

//' Get subset
//'
//' From a list of spectra, get the subset which is described by an index vector.
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @param index an integer vector holding the indices of the spectra to retrieve
//' @return the requested subset of spectra returned as a list
// [[Rcpp::export]]
List getSubset(List L, IntegerVector index) {
    /* the length of our list is the number of spectra */
    int nSpectra = L.length();
    int ndx = index.length();
    // Rcout << "spectra = " << nSpectra << ", length of index = " << ndx << std::endl;
    List out(ndx);

    for (int i = 0; i < ndx; i++) {
        out[i] = L[index[i]];
    }

    return out;
}

//' Calculate average spectrum
//'
//' From a list of spectra, calculate the average spectrum where the weighting is
//' done using system temperature and integration times of the headers.
//'
//' @param L a list of spectra (each with components 'head', 'freq' and 'data')
//' @return the average spectrum
// [[Rcpp::export]]
List average(List L) {
    static char error[80];

    /* get first spectrum */
    List l = L[0];

    Environment env = Environment::global_env();
    List options = env["drp.options"];
    double tol = as<double>(options["position.tolerance"]);
    Rcout << "position tolerance = " << tol << std::endl;

    int nSpectra = L.length();
    /* get 'head' part of that first spectrum */
    List head = l[0];
    /* get 'freq' part of that first spectrum */
    NumericVector freq0 = l[1];
    int nChannels = freq0.length();
    /* get 'data' part of that first spectrum */
    NumericVector data0 = l[2];
    if (freq0.length() != data0.length()) {
        sprintf(error, "frequency and data vectors differ in length: %ld <> %ld",
                freq0.length(), data0.length());
        stop(error);
    }

    double dt = 1.0;
    if (head.containsElementNamed("dt")) dt = as<double>(head["dt"]);

    double Tsys = 300.0;
    if (head.containsElementNamed("T.sys")) Tsys = as<double>(head["T.sys"]);

    double w = dt/(Tsys*Tsys);
    double total = dt*pow(300.0/Tsys, 2.0);

    for (int k = 0; k < nChannels; k++) data0[k] *= w;
    Rcout << "dt = " << dt << ", Tsys = " << Tsys << std::endl;

    double sumw = w;
    for (int i = 1; i < nSpectra; i++) {
        l = L[i];
        head = l[0];
        if (head.containsElementNamed("dt"))    dt = as<double>(head["dt"]);
        if (head.containsElementNamed("T.sys")) Tsys = as<double>(head["T.sys"]);
        w = dt/(Tsys*Tsys);
        sumw += w;
        total += dt*pow(300.0/Tsys, 2.0);

        NumericVector freq = l[1];
        if (freq.length() != freq0.length()) {
            sprintf(error, "frequency vector %d differs in length: %ld <> %ld",
                    i, freq.length(), freq0.length());
            stop(error);
        }

        NumericVector data = l[2];
        if (data.length() != data0.length()) {
            sprintf(error, "data vector %d differs in length: %ld <> %ld",
                    i, data.length(), data0.length());
            stop(error);
        }
        for (int k = 0; k < nChannels; k++) {
            data0[k] += data[k]*w;
        }
    }
    for (int k = 0; k < nChannels; k++) data0[k] /= sumw;

    if (head.containsElementNamed("dt")) head["dt"] = total;
    if (head.containsElementNamed("T.sys")) head["T.sys"] = 300.0;

    List average = List::create(Named("head") = head, Named("freq") = freq0, Named("data") = data0);
    average.attr("class") = "spectrum";

    return average;
}

// [[Rcpp::export]]
void modify(List L, std::string column, SEXP value) {
    static char error[80];

    /* get the first spectrum */
    List l = L[0];
    List head = l[0];

    int nSpectra = L.length();
    int typ = TYPEOF(value);
    int n = Rf_length(value);
    // Rcout << "value of type " << typ << " and length " << n << std::endl;
    if ((n != 1) && (n != nSpectra)) {
        sprintf(error, "value should have length 1 or %d", nSpectra);
        stop(error);
    }

    GenericVector v(value);
    for (int i = 0; i < nSpectra; i++) {
        l = L[i];
        head = l[0];
        if (head.containsElementNamed(column.data())) {
            if (n == 1) head[column] = v[0];
            else        head[column] = v[i];
        }
    }
}

// [[Rcpp::export]]
void addColumns(List L, CharacterVector newnames) {
    static char error[80];

    int nColumns = newnames.length();
    int nSpectra = L.length();
    List l = L[0];
    List head = l[0];
    StringVector names = head.names();
    for (int j = 0; j < nColumns; j++) names.push_back(newnames[j]);

    IntegerVector v(nSpectra);
    for (int i = 0; i < nSpectra; i++) {
        l = L[i];
        head = l[0];
        int oldsize = head.length();
        List newhead(oldsize+nColumns);
        for (int j = 0; j < oldsize; j++) {
            newhead[j] = head[j];
        }
        // for (int j = 0; j < nColumns; j++) {
        //     newhead[j] = v[i];
        // }
        newhead.attr("names") = names;
        l[0] = newhead;
    }
}

// Prototype of function working on a single spectrum
//
// Does nothing
// [[Rcpp::export]]
List foo(List S) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S[0];
    NumericVector freq = S[1];
    NumericVector data = S[2];

    List S1 = List::create(Named("head") = head, Named("freq") = freq, Named("data") = data);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Clip channels from spectra
//'
//' Given a single spectrum and a vector of channel numbers, clip the spectrum such that
//' only the given channels are kept in the frequency and data vectors.
//'
//' @param S a single spectrum
//' @param keep a vector holding the channel numbers to keep
//' @return the clipped spectrum
// [[Rcpp::export]]
List clip(List S, IntegerVector keep) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S[0];
    NumericVector freq0 = S[1];
    NumericVector data0 = S[2];
    int nk = keep.length();
    int nc = data0.length();

    NumericVector freq1(nk);
    NumericVector data1(nk);
    for (int i = 0; i < nk; i++) {
        freq1[i] = freq0[keep[i]];
        data1[i] = data0[keep[i]];
    }

    List S1 = List::create(Named("head") = head, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

/*** R
# S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)), freq=-5:5, data=rnorm(11))
# S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)), freq=-5:5, data=rnorm(11))
# L <- list(S1, S2)
# getHead(list(S1,S2))

files  <- system("ls /home/olberg/R/salsa/FITS/spectrum*.fits", intern = TRUE)
L <- readSALSA(files)

index <- seq(3,11, by=2)
L1 <- getSubset(L, index)
head <- getHead(L1)
F <- getFreq(L1)
D <- getData(L1)
matplot(F, D, lty=1, type='l', xlab="frequency", ylab="")
A <- average(L1)
lines(A$freq, A$data, col="black", lwd=3)
A1 <- clip(A, seq(120,220))
lines(A1$freq, A1$data, col="green", lwd=3)
*/
