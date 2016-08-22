#include <Rcpp.h>

using namespace Rcpp;

//' Get frequency vectors
//'
//' From a list of spectra, get the frequency vectors and return as matrix.
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return a matrix having all the frequency vectors as columns
//' @seealso \code{\link{getVelo}}
// [[Rcpp::export]]
NumericMatrix getFreq(List L) {
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector freq = l["freq"];
    int nChannels = freq.length();

    /* allocate matrix */
    NumericMatrix F(nChannels, nSpectra);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        freq = l["freq"];
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

//' Construct frequency vector
//'
//' Turn a given velocity vector into a frequency vector.
//'
//' @param S a single spectrum
//' @param v a numeric vector holding velocities
//' @return a numeric vector holding frequencies
//' @seealso \code{\link{velocities}}
//' @examples
//' # resample spectrum at particular velocities
//' \dontrun{
//' v <- seq(-100,100,by=2)
//' R <- resample(S, frequencies(S,v))
//' }
//' NULL
// [[Rcpp::export]]
NumericVector frequencies(List S, NumericVector v) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S["head"];

    int nc = v.length();
    NumericVector freq(nc);

    double f0 = 0.0;
    double vs = 0.0;

    /* initialize vs to mean of frequency vector */
    for (int i = 0; i < nc; i++) vs += v[i];
    vs /= nc;

    /* use header variables, when available */
    if (head.containsElementNamed("f0"))    f0 = as<double>(head["f0"]);
    if (head.containsElementNamed("v.LSR")) vs = as<double>(head["v.LSR"]);
    // Rcout << "f0 = " << f0 << ", vs = " << vs << std::endl;

    if (f0 == 0.0) stop("zero frequency");
    /* now do conversion */
    for (int i = 0; i < nc; i++) {
        freq[i] = f0*(1.0-(v[i]-vs)/299792.46);
    }
    return freq;
}

//' Construct velocity vector
//'
//' Turn the frequency vector into a velocity vector.
//'
//' @param S a single spectrum
//' @return a numeric vector holding velocities
//' @seealso \code{\link{frequencies}}
// [[Rcpp::export]]
NumericVector velocities(List S) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S["head"];

    NumericVector freq = S["freq"];
    int nc = freq.length();
    NumericVector velo(nc);

    double f0 = 0.0;
    double vs = 0.0;
    // Rcout << "1 f0 = " << f0 << ", vs = " << vs << std::endl;

    /* initialize f0 to mean of frequency vector */
    for (int i = 0; i < nc; i++) f0 += freq[i];
    f0 /= nc;
    // Rcout << "2 f0 = " << f0 << ", vs = " << vs << std::endl;

    /* use header variables, when available */
    if (head.containsElementNamed("f0"))    f0 = as<double>(head["f0"]);
    // Rcout << "3 f0 = " << f0 << ", vs = " << vs << std::endl;

    if (head.containsElementNamed("v.LSR")) vs = as<double>(head["v.LSR"]);
    // Rcout << "4 f0 = " << f0 << ", vs = " << vs << std::endl;

    if (f0 == 0.0) stop("zero frequency");
    /* now do conversion */
    for (int i = 0; i < nc; i++) {
        velo[i] = (f0-freq[i])*299792.46/f0 + vs;
    }
    return velo;
}

//' Get velocity vectors
//'
//' From a list of spectra, get the velocity vectors and return as matrix.
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return a matrix having all the velcity vectors as columns
//' @seealso \code{\link{getFreq}}
// [[Rcpp::export]]
NumericMatrix getVelo(List L) {
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    NumericMatrix V = getFreq(L);
    int nChannels = V.nrow();

    for (int i = 0; i < nSpectra; i++) {
        List S = L[i];
        NumericVector v = velocities(S);
        for (int j = 0; j < nChannels; j++) {
            V(j, i) = v[j];
        }
    }
    return V;
}

//' Get data vectors
//'
//' From a list of spectra, get the data vectors and return as matrix.
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return a matrix having all the data vectors as columns
//' @seealso \code{\link{getFreq}}, \code{\link{getVelo}}
//' @examples
//' data(salsa)
//' D <- getData(salsa)
//' image(D)             # show data matrix as color image
// [[Rcpp::export]]
NumericMatrix getData(List L) {
    static char error[80];

    /* the length of our list is the number of spectra */
    int nSpectra = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector data = l["data"];
    int nChannels = data.length();

    /* allocate matrix */
    NumericMatrix D(nChannels, nSpectra);

    for (int row = 0; row < nSpectra; row++) {
        /* get spectrum number 'row' */
        l = L[row];
        data = l["data"];
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
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return a two component integer vector (nChannels, nSpectra)
// [[Rcpp::export]]
IntegerVector getDimension(List L) {
    IntegerVector dim(2);

    /* the length of our list is the number of spectra */
    dim[1] = L.length();

    /* get the first spectrum */
    List l = L[0];

    /* get the data part, their length is the number of channels */
    NumericVector data = l["data"];
    dim[0] = data.length();

    return dim;
}

//' Calculate average spectrum
//'
//' From a list of spectra, calculate the average spectrum where the weighting is
//' done using system temperature and integration times of the headers.
//'
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return the average spectrum
//' @examples
//' data(salsa)
//' A <- average(salsa)
//' plot(salsa)                              # plot individual spectra
//' lines(A$freq, A$data, lwd=5, col='red')  # draw average spectrum on top
// [[Rcpp::export]]
List average(List L) {
    static char error[80];

    /* get first spectrum */
    List l = L[0];

    int nSpectra = L.length();
    /* get 'head' part of that first spectrum */
    List head0 = l["head"];
    /* get 'freq' part of that first spectrum */
    NumericVector freq0 = l["freq"];
    int nChannels = freq0.length();
    /* get 'data' part of that first spectrum */
    NumericVector data0 = l["data"];
    if (freq0.length() != data0.length()) {
        sprintf(error, "frequency and data vectors differ in length: %ld <> %ld",
                freq0.length(), data0.length());
        stop(error);
    }

    List head1 = clone(head0);
    NumericVector freq1 = clone(freq0);
    NumericVector data1 = clone(data0);

    double dt = 1.0;
    if (head0.containsElementNamed("dt")) dt = as<double>(head0["dt"]);

    double Tsys = 300.0;
    if (head0.containsElementNamed("T.sys")) Tsys = as<double>(head0["T.sys"]);
    // Rcout << "dt = " << dt << ", Tsys = " << Tsys << std::endl;

    double w = dt/(Tsys*Tsys);
    double total = dt*pow(300.0/Tsys, 2.0);

    for (int k = 0; k < nChannels; k++) data1[k] *= w;

    double sumw = w;
    for (int i = 1; i < nSpectra; i++) {
        l = L[i];
        List head = l["head"];
        if (head.containsElementNamed("dt"))    dt = as<double>(head["dt"]);
        if (head.containsElementNamed("T.sys")) Tsys = as<double>(head["T.sys"]);
        // Rcout << "dt = " << dt << ", Tsys = " << Tsys << std::endl;
        w = dt/(Tsys*Tsys);
        sumw += w;
        total += dt*pow(300.0/Tsys, 2.0);

        NumericVector freq = l["freq"];
        if (freq.length() != freq0.length()) {
            sprintf(error, "frequency vector %d differs in length: %ld <> %ld",
                    i, freq.length(), freq0.length());
            stop(error);
        }

        NumericVector data = l["data"];
        if (data.length() != data1.length()) {
            sprintf(error, "data vector %d differs in length: %ld <> %ld",
                    i, data.length(), data1.length());
            stop(error);
        }
        for (int k = 0; k < nChannels; k++) {
            data1[k] += data[k]*w;
        }
    }
    for (int k = 0; k < nChannels; k++) data1[k] /= sumw;

    if (head1.containsElementNamed("dt"))    head1["dt"] = total;
    if (head1.containsElementNamed("T.sys")) head1["T.sys"] = 300.0;

    List average = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    average.attr("class") = "spectrum";

    return average;
}

void printStats(const char *cmd, NumericVector d) {
    int nc = d.length();
    double min, max, mean, var;

    mean = min = max = d[0];
    var = 0.0;
    for (int i = 1; i < nc; i++) {
        if (min > d[i]) min = d[i];
        if (max < d[i]) max = d[i];
        mean += d[i];
    }
    mean /= (double)nc;
    for (int i = 0; i < nc; i++) var += (d[i]-mean)*(d[i]-mean);
    var /= (double)(nc-1);
    Rcout << cmd << ":" << min << " " << max << " " << mean << " " << var << std::endl;
}

bool useVelocity() {
    bool flag = false;
    static const SEXP option = Rf_install("system");
    std::string system = as<std::string>(Rf_GetOption1(option));

    flag = (system == "velocity");
    return flag;
}

// Prototype of function working on a single spectrum
//
// Does nothing
// [[Rcpp::export]]
List foo(List S) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head0 = S["head"];
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];

    List head1 = clone(head0);
    NumericVector freq1 = clone(freq0);
    NumericVector data1 = clone(data0);

    bool flag = useVelocity();

    // put code here
    
    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Fold a frequency switched spectrum
//'
//' @param S a single spectrum
//' @param ft frequency requency throw in Mhz
//' @param shift if TRUE, assume symmetric switching around true, centre frequency
//' @return the folded spectrum
// [[Rcpp::export]]
List fold(List S, double ft, bool shift = false) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head0 = S["head"];
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];
    int nc = data0.length();

    List head1 = clone(head0);
    NumericVector freq1 = clone(freq0);
    NumericVector data1 = clone(data0);

    if (ft == 0.0) stop("zero frequency throw");
    int n = (int)(ft/(freq0[1]-freq0[0]));
    int i;

    if (n > 0) {
	for (i = 0; i < nc-n; i++) {
	    data1[i] -= data1[i+n];
	    data1[i] /= 2.0;
	}
	if (shift) {
	    n /= 2;
	    for (i = nc-1; i >= n; i--) {
		data1[i] = data1[i-n];
	    }
	}
    } else {
	n = -n;
	for (i = nc-1; i >= n; i--) {
	    data1[i] -= data1[i-n];
	    data1[i] /= 2.0;
	}
	if (shift) {
	    n /= 2;
	    for (i = 0; i < nc-n; i++) {
		data1[i+n] = data1[i];
	    }
	}
    }

    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Reverse a spectrum
//'
//' Reverse the order of the channels, i.e. turn both the frequency and data
//' vector around.
//' @param S a single spectrum
//' @return the reversed spectrum
// [[Rcpp::export]]
List reverse(List S) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head0 = S["head"];
    List head1 = clone(head0);
    if (head1.containsElementNamed("df")) {
        double df = as<double>(head1["df"]);
        head1["df"] = -df;
    }
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];

    int nc = data0.length();
    NumericVector freq1(nc);
    NumericVector data1(nc);
    // Rprintf("old: %p %p %p\n", head0, freq0, data0);
    // Rprintf("new: %p %p %p\n", head1, freq1, data1);

    for (int i = 0; i < nc; i++) {
        freq1[i] = freq0[nc-1-i];
        data1[i] = data0[nc-1-i];
    }
    // printStats("reverse", data1);

    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Calculate integrated area
//'
//' Given a mask defining the spectral areas of interest, return the integrated
//' area over those region(s).
//' @param S a single spectrum
//' @param mask a logical vector equal to TRUE for all the channels that should
//'        be integrated.
//' @return the integrated value
//' @examples
//' data(salsa)
//' options(system="velocity")    # work in velocity space
//' S <- salsa[[1]]               # get the first spectrum
//' v <- velocities(S)            # get velocity vector
//' mask <- (v > -20) & (v < 20)  # integrate from -20..20 km/s
//' area(S, mask)                 # calculate integrated area in K*km/s
//' # call 'area' for each of the spectra in 'salsa' with parameter 'mask'
//' sapply(salsa, FUN=area, mask)
// [[Rcpp::export]]
double area(List S, LogicalVector mask) {
    static char error[80];

    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S["head"];
    NumericVector freq = S["freq"];
    NumericVector data = S["data"];
    int nc = data.length();
    if (mask.length() != nc) {
        sprintf(error, "supplied mask has wrong length: %ld <> %d", mask.length(), nc);
        stop(error);
    }

    bool flag = useVelocity();
    double dx = 0.0;
    double sum = 0.0;
    for (int i = 1; i < nc-1; i++) {
        if (mask[i]) {
            dx = fabs(freq[i+1]-freq[i-1])/2.0;
            if (flag) dx = 2.0*299792.46*dx/(freq[i+1]+freq[i-1]);
            sum += data[i]*dx;
        }
    }

    return sum;
}

//' Calculate moments
//'
//' Given a mask defining the spectral areas of interest, return the integrated
//' area over those region(s).
//' @param S a single spectrum
//' @param mask a logical vector equal to TRUE for all the channels that should
//'        be integrated.
//' @return a data frame containing minimum, maximum, mean, sdev, skewness and kurtosis
//' @seealso \code{\link{area}}
//' @examples
//' data(salsa)
//' A <- average(salsa)
//' lmask <- mask(A, c(1420,1421))  # calculate moments between 1420...1421 MHz
//' moment(A, mask=lmask)
//' do.call("rbind", lapply(salsa, moment, lmask)) # do this for individual spectra
// [[Rcpp::export]]
DataFrame moment(List S, LogicalVector mask) {
    static char error[80];

    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head = S["head"];
    NumericVector freq = S["freq"];
    NumericVector data = S["data"];
    int nc = data.length();
    if (mask.length() != nc) {
        sprintf(error, "supplied mask has wrong length: %ld <> %d", mask.length(), nc);
        stop(error);
    }

    double mean = data[0];
    double Tmin = mean;
    double Tmax = Tmin;
    double var = 0.0;
    double sdev = 0.0;
    double skew = 0.0;
    double kurt = 0.0;
    for (int i = 1; i < nc; i++) {
        double X = data[i];
        if (X < Tmin) Tmin = X;
        if (X > Tmax) Tmax = X;
        mean += X;
    }
    mean /= (double)nc;

    for (int i = 0; i < nc; i++) {
        double X = data[i];
        double Y = X*X;
        var += Y;
        Y *= X;
        skew += Y;
        Y *= X;
        kurt += Y;
    }
    var /= (double)(nc-1);
    sdev = sqrt(var);
    if (var != 0.0) {
        skew /= pow((double)nc*sdev, 3.0);
        kurt /= pow((double)nc*var, 2.0)-3.0;

    }
    return DataFrame::create(_["min"]= Tmin, _["max"]= Tmax, _["mean"]= mean,
                             _["sdev"]= sdev, _["skew"]= skew, _["kurt"]= kurt);
}

//' Trim channels from spectra
//'
//' Given a single spectrum and a vector of channel numbers, trim the spectrum such that
//' only the given channels are kept in the frequency and data vectors.
//'
//' @param S a single spectrum
//' @param keep a vector holding the channel numbers to keep
//' @return the trimmed spectrum
// [[Rcpp::export]]
List trim(List S, IntegerVector keep) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");
    List head0 = S["head"];
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];
    int nk = keep.length();
    int nc = data0.length();

    List head1 = clone(head0);
    NumericVector freq1(nk);
    NumericVector data1(nk);

    // Rprintf("old: %p %p %p\n", head0, freq0, data0);
    // Rprintf("new: %p %p %p\n", head1, freq1, data1);

    for (int i = 0; i < nk; i++) {
        freq1[i] = freq0[keep[i]-1];
        data1[i] = data0[keep[i]-1];
    }
    // printStats("trim", data1);

    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Filter spectrum
//'
//' Perform filtering on a single spectrum.
//'
//' @param S a single spectrum
//' @param coeffs a numeric vector with an odd number of filter coefficients
//' @return the filtered spectrum
// [[Rcpp::export]]
List sieve(List S, NumericVector coeffs) {
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");

    List head0 = S["head"];
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];
    int nf = coeffs.length();
    if ((nf % 2) == 0) stop("number of coefficients must be odd");

    int nc = data0.length();

    /* make sure filter coefficients are normalized */
    double sum = 0.0;
    for (int j = 0; j < nf; j++) sum += coeffs[j];
    for (int j = 0; j < nf; j++) coeffs[j] /= sum;

    /* calculate length of output vectors */
    List head1 = clone(head0);
    nf = (nf+1)/2;
    int nout = (nc-2*nf)/nf + 1;
    NumericVector freq1(nout);
    NumericVector data1(nout);

    // Rprintf("old: %p %p %p\n", head0, freq0, data0);
    // Rprintf("new: %p %p %p\n", head1, freq1, data1);

    int iout = 0;
    for (int i = nf-1; i < nc-nf; i += nf) {
        freq1[iout] = freq0[i];
        data1[iout] = 0;
        for (int j = -(nf-1); j <= (nf-1); j++) {
            data1[iout] += data0[i+j]*coeffs[j+nf-1];
            // Rcout << "data0[" << i+j << "] = " << data0[i+j] << std::endl;
            // Rcout << "coeffs[" << j+nf-1 << "] = " << coeffs[j+nf-1] << std::endl;
            // Rcout << "data1[" << iout << "] = " << data1[iout] << std::endl;
        }
        iout++;
    }
    // printStats("sieve", data1);

    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Resample spectrum
//'
//' Perform resampling (cubic spline interpolation) on a single spectrum.
//'
//' @param S a single spectrum
//' @param f a frequency vector onto which the spectrum should be resampled
//' @param smooth if TRUE convolve with Gaussian response
//' @return the resampled spectrum
// [[Rcpp::export]]
List resample(List S, NumericVector f, bool smooth=false) {
    // We will use a Gaussian filter below with a noise bandwidth given by the
    // channel spacing of vector f, we need the follwoing constants for this:
    const double a = 4.0*log(2.0);
    const double A = sqrt(M_PI/a);
    double df, DF, Df, df0, df1, w, sum, sumw;
    int i, j;
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");

    List head0 = S["head"];
    List head1 = clone(head0);
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];

    int nc = data0.length();
    int nout = f.length();
    // Rcout << "resample: f[0] = " << f[0] << ", f[" << nout-1 << "] = " << f[nout-1] << std::endl;
    NumericVector freq1 = clone(f);
    NumericVector data1(nout);

    df = fabs(freq0[1]-freq0[0]);
    DF = fabs(freq1[1]-freq1[0]);
    if (DF < df) stop("Can't resample to higher resolution");
    df0 = (DF-df)/2.0;
    df1 = (DF+df)/2.0;

    if (smooth) {
        for (j = 0; j < nout; j++) {
            if ((freq1[j] < freq0[0]) || (freq1[j] > freq0[nc-1])) {
                data1[j] = NA_REAL;
            } else {
                sum = 0.0;
                sumw = 0.0;
                w = 0.0;
                for (i = 0; i < nc; i++) {
                    double delta = A*(freq1[j]-freq0[i])/DF;
                    if (fabs(delta) < 1.5) {
                        w = exp(-a*delta*delta);
                        sum += data0[i]*w;
                        sumw += w;
                    }
                }
                data1[j] = sum/sumw;
            }
        }
    } else {
        for (j = 0; j < nout; j++) {
            if ((freq1[j] < freq0[0]) || (freq1[j] > freq0[nc-1])) {
                data1[j] = NA_REAL;
            } else {
                sum = 0.0;
                sumw = 0.0;
                w = 0.0;
                for (i = 0; i < nc; i++) {
                    Df = fabs(freq1[j] - freq0[i]);
                    if (Df < df1) {
                        if (Df < df0) w = 1.0;
                        else          w = (df1-Df)/(df1-df0);
                        sum += data0[i]*w;
                        sumw += w;
                    }
                }
                data1[j] = sum/sumw;
            }
        }
    }

    if (head1.containsElementNamed("df"))    head1["df"] = DF;
    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Rescale spectrum
//'
//' Add bias and/or scale by factor.
//'
//' @param S a single spectrum
//' @param factor a numeric value by which to scale the whole spectrum
//' @param bias a numeric value to add to all channels
//' @return the rescaled spectrum, out = in * factor + bias
// [[Rcpp::export]]
List rescale(List S, double factor = 1.0, double bias = 0.0) {
    int i, j;
    if (!S.inherits("spectrum")) stop("Input must be a spectrum");

    List head0 = S["head"];
    NumericVector freq0 = S["freq"];
    NumericVector data0 = S["data"];

    List head1 = clone(head0);
    NumericVector freq1 = clone(freq0);
    NumericVector data1 = clone(data0);

    int nc = data0.length();

    for (i = 1; i < nc; i++) {
        data1[i] = data0[i]*factor + bias;
    }

    List S1 = List::create(Named("head") = head1, Named("freq") = freq1, Named("data") = data1);
    S1.attr("class") = "spectrum";
    return S1;
}

//' Construct a line mask
//'
//' For a given set of frequency (or velocity) windows, construct a
//' matrix of logical values, defining the areas which should not be used
//' in baseline fitting.
//' @param S a single spectrum
//' @param limits pairs of values, which each define a window
//' @return a vector of logical values, one per channel
// [[Rcpp::export]]
LogicalVector mask(List S, NumericVector limits) {
    int nw = limits.length();
    if ((nw % 2) == 1) {
        stop("length of limits must be even");
    }
    nw /= 2;
    NumericVector freq = S["freq"];
    int nc = freq.length();

    bool flag = useVelocity();
    NumericVector x(clone(freq));
    if (flag) x = velocities(S);

    LogicalVector window(nc);
    for (int j = 0; j < nc; j++) {
        window[j] = false;
        for (int i = 0; i < nw; i++) {
            if ((x[j] >= limits[2*i]) && (x[j] <= limits[2*i+1])) {
                window[j] = true;
                break;
            }
        }
    }
    return window;
}

/*** R
# S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)),
#            freq=-5:5, data=rnorm(11))
# S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)),
#            freq=-5:5, data=rnorm(11))
# L <- list(S1, S2)
# getHead(L)

files  <- system("ls /home/olberg/R/salsa/FITS/spectrum*.fits", intern = TRUE)
L <- readSALSA(files)

index <- seq(3,11, by=2)
L1 <- L[index]
head <- getHead(L1)
F <- getFreq(L1)
D <- getData(L1)
matplot(F, D, lty=1, type='l', xlab="frequency", ylab="")
A <- average(L1)
lines(A$freq, A$data, col="black", lwd=3)

M <- mask(A, range(A$freq[120:220]))
area(A, M)
lines(A$freq, M*10, type='S', lwd=3, col='grey')
A1 <- trim(A, seq(120,220))
lines(A1$freq, A1$data, col="green", lwd=3)
*/
