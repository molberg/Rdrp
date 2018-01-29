#include <R.h>
#include <Rinternals.h>

//' Get header data.
//'
//' Given a list L, where each list member is itself a list with
//' components 'head' (which is a list or a dataframe with one row),
//' 'freq' (a numeric vector) and 'data' (another numeric vector of same
//' length as freq), return all the 'head' components as a data frame,
//' which will have as many rows as the length of the original list.
//'
//' It is assumed that all 'head' components have the same number of
//' members with identical names and types.
//'
//' @param L a list of spectra, each with components 'head', 'freq' and 'data'
//' @return a data.frame formed by row-binding all the individual 'head's
//' @seealso \code{\link{modify}}
//' @examples
//' S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)),
//'            freq=-5:5, data=rnorm(11))
//' S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)),
//'            freq=-5:5, data=rnorm(11))
//' L <- list(S1,S2)
//' class(L) <- "spectra"
//' getHead(L)
//'
//' # will result in
//'
//' #   target   ra  dec dt
//' # 1  Orion 1.23 -0.5 20
//' # 2  SgrB2 5.43  0.5 20
//'
//' @export
// [[Rcpp::export]]
SEXP getHead(SEXP L) {
    SEXP col;
    int nx = Rf_length(L);
    SEXP cls = Rf_getAttrib(L, R_ClassSymbol);
    int isSpectra = strncmp(CHAR(asChar(cls)), "spectra", 7);
    if (isSpectra != 0) {
        warning("list is not of class 'spectra'\n");
        return R_NilValue;
    }
    // printf("L is of class: %s\n", CHAR(asChar(cls)));
    // Rprintf("nx = %d\n", nx);
    SEXP S = VECTOR_ELT(L, 0);
    SEXP head = VECTOR_ELT(S, 0);
    int nh = Rf_length(head);
    // Rprintf("nh = %d\n", nh);
    SEXP H = PROTECT(allocVector(VECSXP, nh));
    SEXP row = PROTECT(allocVector(INTSXP, nx));
    namesgets(H, getAttrib(head, R_NamesSymbol));
    for (int i = 0; i < nh; i++) {
        INTEGER(row)[i] = i+1;
        SEXP col = VECTOR_ELT(head, i);
        int typ = TYPEOF(col);
        // Rprintf("col %d typ = %d\n", i+1, typ);
        switch(typ) {
          case INTSXP: // integer
            col = PROTECT(allocVector(INTSXP, nx));
            for (int j = 0; j < nx; j++) {
                S = VECTOR_ELT(L, j);
                head = VECTOR_ELT(S, 0);
                INTEGER(col)[j] = asInteger(VECTOR_ELT(head, i));
            }
            SET_VECTOR_ELT(H, i, col);
            UNPROTECT(1);
            break;
          case REALSXP: // float
            col = PROTECT(allocVector(REALSXP, nx));
            for (int j = 0; j < nx; j++) {
                S = VECTOR_ELT(L, j);
                head = VECTOR_ELT(S, 0);
                REAL(col)[j] = asReal(VECTOR_ELT(head, i));
            }
            SET_VECTOR_ELT(H, i, col);
            UNPROTECT(1);
            break;
          case STRSXP: // string
            col = PROTECT(allocVector(STRSXP, nx));
            for (int j = 0; j < nx; j++) {
                S = VECTOR_ELT(L, j);
                head = VECTOR_ELT(S, 0);
                SET_STRING_ELT(col, j, mkChar(CHAR(asChar(VECTOR_ELT(head, i)))));
            }
            SET_VECTOR_ELT(H, i, col);
            UNPROTECT(1);
            break;
          default:
            warning("unsupported column type %d\\n", typ);
            break;
        }
    }
    Rf_setAttrib(H, R_RowNamesSymbol, row);
    UNPROTECT(1);
    Rf_setAttrib(H, R_ClassSymbol, Rf_mkString("data.frame"));
    UNPROTECT(1);
    return H;
}
