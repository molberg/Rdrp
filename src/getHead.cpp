#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>

//' Get one column of header data.
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
//' @param index the index of the column to be retrieved
//' @return a data.frame formed by row-binding all the individual 'head's
//' @seealso \code{\link{modify}}
//' @examples
//' S1 <- list(head=list(target="Orion", ra=1.23, dec=-0.5, dt=as.integer(20)),
//'            freq=-5:5, data=rnorm(11))
//' S2 <- list(head=list(target="SgrB2", ra=5.43, dec=+0.5, dt=as.integer(20)),
//'            freq=-5:5, data=rnorm(11))
//' L <- list(S1,S2)
//' class(L) <- "spectra"
//' i <- which(names(S1$head) == "ra")
//' getHeaderColumn(L, i)
//'
//' # will result in
//'
//' # [1] 1.23 5.43
//'
//' @export
// [[Rcpp::export]]
SEXP getHeaderColumn(SEXP L, SEXP index) {
    SEXP col;
    int idx = INTEGER_VALUE(index);
    const char *cls = STRING_VALUE(Rf_getAttrib(L, R_ClassSymbol));
    int isSpectra = strncmp(cls, "spectra", 7);
    if (isSpectra != 0) {
        Rf_warning("list is not of class 'spectra'\n");
        return R_NilValue;
    }

    /* get length of input list, i.e number of spectra */
    int nx = Rf_length(L);

    /* get first spectrum */
    SEXP S = VECTOR_ELT(L, 0);
    /* get header of first spectrum */
    SEXP head = VECTOR_ELT(S, 0);
    /* get length of header, i.e. number of list elements */
    int nh = Rf_length(head);
    if ((idx < 1) || (idx > nh)) {
        Rf_warning("column index out of range\n");
        return R_NilValue;
    }

    int i = idx-1;
    /* get header element i */
    SEXP next = VECTOR_ELT(head, i);
    /* query its type */
    int typ = TYPEOF(next);

    /* now, based on type construct vectors for each column */
    switch(typ) {
      case LGLSXP: // logical
        col = PROTECT(Rf_allocVector(LGLSXP, nx));
        for (int j = 0; j < nx; j++) {
            SEXP s = VECTOR_ELT(L, j);
            SEXP h = VECTOR_ELT(s, 0);
            SEXP item = VECTOR_ELT(h, i);
            LOGICAL(col)[j] = LOGICAL_VALUE(item);
            /* LOGICAL(col)[j] = 0; */
        }
        break;
      case INTSXP: // integer
        col = PROTECT(Rf_allocVector(INTSXP, nx));
        for (int j = 0; j < nx; j++) {
            SEXP s = VECTOR_ELT(L, j);
            SEXP h = VECTOR_ELT(s, 0);
            SEXP item = VECTOR_ELT(h, i);
            INTEGER(col)[j] = INTEGER_VALUE(item);
            /* INTEGER(col)[j] = j; */
        }
        break;
      case REALSXP: // float
        col = PROTECT(Rf_allocVector(REALSXP, nx));
        for (int j = 0; j < nx; j++) {
            SEXP s = VECTOR_ELT(L, j);
            SEXP h = VECTOR_ELT(s, 0);
            SEXP item = VECTOR_ELT(h, i);
            REAL(col)[j] = NUMERIC_VALUE(item);
            /* REAL(col)[j] = 2.0*j; */
        }
        break;
      case STRSXP: // string
        col = PROTECT(Rf_allocVector(STRSXP, nx));
        for (int j = 0; j < nx; j++) {
            SEXP s = VECTOR_ELT(L, j);
            SEXP h = VECTOR_ELT(s, 0);
            SEXP item = VECTOR_ELT(h, i);
            SET_STRING_ELT(col, j, Rf_mkChar(STRING_VALUE(item)));
            /* SET_STRING_ELT(col, j, mkChar("NA")); */
        }
        break;
      default:
        Rf_warning("unsupported column type %d\\n", typ);
        break;
    }
    UNPROTECT(1); // col
    return col;
}

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

    /* check that input list is of class 'spectra' */
    const char *cls = STRING_VALUE(Rf_getAttrib(L, R_ClassSymbol));
    int isSpectra = strncmp(cls, "spectra", 7);
    if (isSpectra != 0) {
        Rf_warning("list is not of class 'spectra'\n");
        return R_NilValue;
    }

    /* get length of input list, i.e number of spectra */
    int nx = Rf_length(L);

    /* get first spectrum */
    SEXP S = VECTOR_ELT(L, 0);
    /* get header of first spectrum */
    SEXP head = VECTOR_ELT(S, 0);
    /* get length of header, i.e. number of list elements */
    int nh = Rf_length(head);
    // Rprintf("nx * nh = %d * %d\n", nx, nh);

    /* copy the column names from the first header */
    SEXP H = PROTECT(Rf_allocVector(VECSXP, nh));
    Rf_namesgets(H, Rf_getAttrib(head, R_NamesSymbol));

    for (int i = 0; i < nh; i++) {
        /* use the row count as the name of the next row */
        /* INTEGER(row)[i] = i+1; */
        col = PROTECT(getHeaderColumn(L, Rf_ScalarInteger(i+1)));
        SET_VECTOR_ELT(H, i, col);
    }
    /* release all allocated column vectors */
    UNPROTECT(nh);

    /*
    for (int i = 0; i < nh; i++) {
        int typ = TYPEOF(VECTOR_ELT(H, i));
        int len = Rf_length(VECTOR_ELT(H, i));
        Rprintf("type of column %d: %d, length %d\n", i, typ, len);
        Rf_PrintValue(VECTOR_ELT(H, i));
    }
    */

    SEXP rnms = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(rnms)[0] = NA_INTEGER;
    INTEGER(rnms)[1] = -nx;

    /* set row names and class attributes of our data frame */
    // Rprintf("setting attributes\n");
    Rf_setAttrib(H, R_ClassSymbol, Rf_ScalarString(Rf_mkChar("data.frame")));
    Rf_setAttrib(H, R_RowNamesSymbol, rnms);

    /* Rf_setAttrib(H, R_ClassSymbol, Rf_mkString("data.frame")); */
    /* Rf_setAttrib(H, R_RowNamesSymbol, row); */
    UNPROTECT(2); // rnms, H
    return H;
}
