#include "class.h"

#undef DEBUG

typedef struct {
    int origpos;
    const char *value;
} SORT;

int fileType(const char *cfname)
{
    FILE *cfp = fopen(cfname, "r");
    if (!cfp) return -1;

    char code[4];
    size_t len = fread(code, sizeof(char), 4, cfp);
    if (len != 4) return -2;

    code[2] = '\0';
    code[3] = '\0';
    fclose(cfp);

    if (strncmp((const char *)code, "1A", 2) == 0) return 1;
    if (strncmp((const char *)code, "2A", 2) == 0) return 2;
    return 0;
}

int qcmp(const void *x, const void *y) {
    int res = strcmp(((SORT *)x)->value, ((SORT *)y)->value);
    if (res != 0)
        return res;
    else
        // they are equal - use original position as tie breaker
        return (((SORT *)x)->origpos - ((SORT *)y)->origpos);
}

SEXP makeFactor(SEXP str)
{
    SORT *sorted;
    int i, level;
    int num;
    const char *ptr;

    num = Rf_length(str);
    sorted = (SORT *)malloc(num*sizeof(SORT));

    for (i = 0; i < num; i++) {
        ptr = CHAR(STRING_ELT(str, i));
        sorted[i].value = ptr;
        sorted[i].origpos = i;
    }
    qsort(sorted, num, sizeof(SORT), qcmp);

    // find number of levels
    level = 1;
    for (i = 0; i < num - 1; i++) {
        if (strcmp(sorted[i].value, sorted[i+1].value) != 0) {
            level++;
        }
    }

    SEXP factor = PROTECT(allocVector(INTSXP, num));
    SEXP unique = PROTECT(allocVector(STRSXP, level));
    level = 0;
    INTEGER(factor)[0] = level+1;
    SET_STRING_ELT(unique, 0, mkChar(sorted[0].value));
    for (i = 0; i < num-1; i++) {
        if (strcmp(sorted[i].value, sorted[i+1].value) != 0) {
            level++;
            SET_STRING_ELT(unique, level, mkChar(sorted[i+1].value));
        }
        INTEGER(factor)[i+1] = level+1;
    }
    Rf_setAttrib(factor, R_LevelsSymbol, unique);
    Rf_setAttrib(factor, R_ClassSymbol, Rf_mkString("factor"));
    free(sorted);

    UNPROTECT(2); // factor, unique
    return factor;
}

void makePOSIXct(SEXP &t)
{
    SEXP cl = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(cl, 0, Rf_mkChar("POSIXct"));
    SET_STRING_ELT(cl, 1, Rf_mkChar("POSIXt"));
    Rf_setAttrib(t, R_ClassSymbol, cl);
    UNPROTECT(1); // cl

    SEXP tz = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(tz, 0, mkChar("UTC"));
    setAttrib(t, Rf_install("tzone"), tz);
    UNPROTECT(1); // tz
}

SEXP emptyFrame(int nspec)
{
    SEXP frame = PROTECT(allocVector(VECSXP, 13));
    SEXP nam = PROTECT(allocVector(STRSXP, 13)); // names attribute (column names)
    SET_STRING_ELT(nam, 0, mkChar("id"));
    SET_STRING_ELT(nam, 1, mkChar("scan"));
    SET_STRING_ELT(nam, 2, mkChar("target"));
    SET_STRING_ELT(nam, 3, mkChar("line"));
    SET_STRING_ELT(nam, 4, mkChar("RA"));
    SET_STRING_ELT(nam, 5, mkChar("Dec"));
    SET_STRING_ELT(nam, 6, mkChar("f.LO"));
    SET_STRING_ELT(nam, 7, mkChar("f0"));
    SET_STRING_ELT(nam, 8, mkChar("df"));
    SET_STRING_ELT(nam, 9, mkChar("v.LSR"));
    SET_STRING_ELT(nam,10, mkChar("dt"));
    SET_STRING_ELT(nam,11, mkChar("T.sys"));
    SET_STRING_ELT(nam,12, mkChar("observed.date"));
    namesgets(frame, nam);
    UNPROTECT(1); // nam

    SEXP id = PROTECT(allocVector(INTSXP, nspec));         // 0
    SEXP scanno = PROTECT(allocVector(INTSXP, nspec));     // 1
    SEXP target = PROTECT(allocVector(STRSXP, nspec));     // 2
    SEXP line = PROTECT(allocVector(STRSXP, nspec));       // 3
    SEXP RA = PROTECT(allocVector(REALSXP, nspec));        // 4
    SEXP Dec = PROTECT(allocVector(REALSXP, nspec));       // 5
    SEXP fLO = PROTECT(allocVector(REALSXP, nspec));       // 6
    SEXP f0 = PROTECT(allocVector(REALSXP, nspec));        // 7
    SEXP df = PROTECT(allocVector(REALSXP, nspec));        // 8
    SEXP vs = PROTECT(allocVector(REALSXP, nspec));        // 9
    SEXP dt = PROTECT(allocVector(REALSXP, nspec));        // 10
    SEXP tsys = PROTECT(allocVector(REALSXP, nspec));      // 11
    SEXP utc = PROTECT(allocVector(REALSXP, nspec));        // 12
    makePOSIXct(utc);
    Rf_setAttrib(frame, R_RowNamesSymbol, id);

    SET_VECTOR_ELT(frame, 0, id);
    SET_VECTOR_ELT(frame, 1, scanno);
    SET_VECTOR_ELT(frame, 2, target);
    SET_VECTOR_ELT(frame, 3, line);
    SET_VECTOR_ELT(frame, 4, RA);
    SET_VECTOR_ELT(frame, 5, Dec);
    SET_VECTOR_ELT(frame, 6, fLO);
    SET_VECTOR_ELT(frame, 7, f0);
    SET_VECTOR_ELT(frame, 8, df);
    SET_VECTOR_ELT(frame, 9, vs);
    SET_VECTOR_ELT(frame,10, dt);
    SET_VECTOR_ELT(frame,11, tsys);
    SET_VECTOR_ELT(frame,12, utc);
    UNPROTECT(13); // id, scanno, target, line, RA, Dec, fLO, f0, df, vs, dt, tsys, utc
    Rf_setAttrib(frame, R_ClassSymbol, Rf_mkString("data.frame"));
    UNPROTECT(1);  // frame

    return frame;
}

SEXP makeSpectrum(SEXP head, SEXP freq, SEXP data)
{
    SEXP S = PROTECT(allocVector(VECSXP, 3)); // head, freq, data

    SEXP nam = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nam, 0, mkChar("head"));
    SET_STRING_ELT(nam, 1, mkChar("freq"));
    SET_STRING_ELT(nam, 2, mkChar("data"));
    namesgets(S, nam);
    UNPROTECT(1); // nam

    setAttrib(S, R_ClassSymbol, Rf_mkString("spectrum"));

    SET_VECTOR_ELT(S, 0, head);
    SET_VECTOR_ELT(S, 1, freq);
    SET_VECTOR_ELT(S, 2, data);
    UNPROTECT(1); // S

    return S;
}

//' Get header information from a GILDAS/CLASS single dish data file
//'
//' Given a filename, open the file and scan it for single dish spectra
//' or continuum scans. A data frame is returned where each row corresponds
//' to the header information of one scan in the file.
//' @param filename name of the GILDAS/CLASS file including path to be opened
//' @return data frame with n rows, where n is the number of scans found
//' @examples
//' filename = "mydata.apex"
//' H <- getClassHeader(filename)
//' print(H)           # print the resulting data frame
//' print(summary(H))  # print a summary
//' @export
// [[Rcpp::export]]
SEXP getClassHeader(SEXP filename)
{
    ClassReader *reader = 0;

    const char *cfname = CHAR(STRING_ELT(filename, 0));

    int type = fileType(cfname);
    if (type == -1) {
        warning("failed to open file '%s'\n", cfname);
        return R_NilValue;
    }
    if (type == -2) {
        warning("failed to determine file type of '%s'\n", cfname);
        return R_NilValue;
    }

    if ((type != 1) && (type != 2)) {
        warning("unrecognized file type!\n");
        return R_NilValue;
    }
    if (type == 1) reader = new Type1Reader(cfname);
    if (type == 2) reader = new Type2Reader(cfname);

    reader->open();
    reader->getFileDescriptor();

    int nscans = reader->getDirectory();
#ifdef DEBUG
    Rprintf("%s: number of scans: %d\n", __FUNCTION__, nscans);
#endif

    SEXP frame = PROTECT(emptyFrame(nscans));
    SEXP id = PROTECT(VECTOR_ELT(frame, 0));
    SEXP scanno = PROTECT(VECTOR_ELT(frame, 1));
    SEXP target = PROTECT(VECTOR_ELT(frame, 2));
    SEXP line = PROTECT(VECTOR_ELT(frame, 3));
    SEXP RA = PROTECT(VECTOR_ELT(frame, 4));
    SEXP Dec = PROTECT(VECTOR_ELT(frame, 5));
    SEXP fLO = PROTECT(VECTOR_ELT(frame, 6));
    SEXP f0 = PROTECT(VECTOR_ELT(frame, 7));
    SEXP df = PROTECT(VECTOR_ELT(frame, 8));
    SEXP vs = PROTECT(VECTOR_ELT(frame, 9));
    SEXP dt = PROTECT(VECTOR_ELT(frame, 10));
    SEXP tsys = PROTECT(VECTOR_ELT(frame, 11));
    SEXP utc = PROTECT(VECTOR_ELT(frame, 12));

    for (int iscan = 0; iscan < nscans; iscan++) {
        SEXP S = PROTECT(reader->getSpectrum(iscan+1, true));

        INTEGER(id)[iscan] = INTEGER(VECTOR_ELT(S, 0))[0];
        INTEGER(scanno)[iscan] = INTEGER(VECTOR_ELT(S, 1))[0];
        SET_STRING_ELT(target, iscan, mkChar(CHAR(STRING_ELT(VECTOR_ELT(S, 2), 0))));
        SET_STRING_ELT(line, iscan , mkChar(CHAR(STRING_ELT(VECTOR_ELT(S, 3), 0))));
        REAL(RA)[iscan] = REAL(VECTOR_ELT(S, 4))[0];
        REAL(Dec)[iscan] = REAL(VECTOR_ELT(S, 5))[0];
        REAL(fLO)[iscan] = REAL(VECTOR_ELT(S, 6))[0];
        REAL(f0)[iscan] = REAL(VECTOR_ELT(S, 7))[0];
        REAL(df)[iscan] = REAL(VECTOR_ELT(S, 8))[0];
        REAL(vs)[iscan] = REAL(VECTOR_ELT(S, 9))[0];
        REAL(dt)[iscan] = REAL(VECTOR_ELT(S, 10))[0];
        REAL(tsys)[iscan] = REAL(VECTOR_ELT(S, 11))[0];
        REAL(utc)[iscan] = REAL(VECTOR_ELT(S, 12))[0];

        UNPROTECT(1); // S
    }
    UNPROTECT(13); // id, scanno, target, line, RA, Dec, f0, fLO, df, vs, dt, tsys, utc;

    SEXP column = PROTECT(VECTOR_ELT(frame, 2));
    target = PROTECT(makeFactor(column));
    SET_VECTOR_ELT(frame, 2, target);
    UNPROTECT(2); // target, column

    column = PROTECT(VECTOR_ELT(frame, 3));
    line = PROTECT(makeFactor(column));
    SET_VECTOR_ELT(frame, 3, line);
    UNPROTECT(2); // line, column

    UNPROTECT(1); // frame
    delete reader;

    return frame;
}

//' Read a GILDAS/CLASS single dish data file
//'
//' Given a filename of a GILDAS/CLASS data file and a data frame 'header',
//' open the file and scan it for single dish spectra or continuum scans.
//' For each row in the header frame, the corresponding scan will be returned
//' as a list with a head, freq and data section. For continuum scans, the
//' frequency vector will simply be a running index. All the individual lists are
//' combined into one major list, which is returned. If you supplied a header
//' and run getHead(L) on the returned list L, you should get your header back.
//' If a header is not supplied, it will be constructed on the fly, such that
//' all scans present in the CLASS file will be returned. Note that this may
//' result in running out of memory for very large (several Gb) files.
//' @param filename name of the GILDAS/CLASS file including path to be opened
//' @param H a data frame with one row for each scan requested (optional).
//' @return list of length n, where n is the number of scans.
//' @examples
//' \dontrun{
//' filename = "mydata.apex"
//' H <- getClassHeader(filename)
//' L = readClass(filename)      # return all spectra in file
//' H <- H[which(H$target == "RDor" & H$line == "CO (6-5)"),]
//' L = readClass(filename, H)   # return only RDor CO (6-5) spectra
//' }
//' @export
// [[Rcpp::export]]
SEXP readClass(SEXP filename, SEXP H = R_NilValue)
{
    ClassReader *reader = 0;

    const char *cfname = CHAR(STRING_ELT(filename, 0));

    int type = fileType(cfname);
    if (type == -1) {
        warning("failed to open file '%s'\n", cfname);
        return R_NilValue;
    }
    if (type == -2) {
        warning("failed to determine file type of '%s'\n", cfname);
        return R_NilValue;
    }

    if ((type != 1) && (type != 2)) {
        warning("unrecognized file type!\n");
        return R_NilValue;
    }
    if (type == 1) reader = new Type1Reader(cfname);
    if (type == 2) reader = new Type2Reader(cfname);

    reader->open();
    reader->getFileDescriptor();

    SEXP id;
    bool all = isNull(H);
    int *index;
    if (all) {
        Rprintf("retrieving all spectra ...\n");
        int nscans = reader->getDirectory();
#ifdef DEBUG
        Rprintf("%s: number of scans: %d\n", __FUNCTION__, nscans);
#endif
        id = PROTECT(allocVector(INTSXP, nscans));
        index = INTEGER(id);
        for (int i = 0; i < nscans; i++) {
            index[i] = i+1;
        }
    } else {
        id = PROTECT(VECTOR_ELT(H, 0));
        index = INTEGER(id);
    }

    int nrows = Rf_length(id);
#ifdef DEBUG
    Rprintf("%s: look-up %d header rows, first = %d, last = %d\n",
            __FUNCTION__, nrows, index[0], index[nrows-1]);
#endif
    SEXP spectra = PROTECT(allocVector(VECSXP, nrows));
    setAttrib(spectra, R_ClassSymbol, Rf_mkString("spectra"));

    for (int i = 0; i < nrows; i++) {
        SEXP S = PROTECT(reader->getSpectrum(index[i], false));
        SET_VECTOR_ELT(spectra, i, S);
        UNPROTECT(1); // S
    }

    UNPROTECT(1); // spectra
    UNPROTECT(1); // id
    delete reader;
    return spectra;
}

ClassReader::ClassReader(const char *filename)
{
    strncpy(cfname, filename, 255);
    m_reclen = 0;
    m_nspec = 0;
}

ClassReader::~ClassReader()
{
    if (cfp != NULL) {
        fclose(cfp);
        cfp = NULL;
    }
}

bool ClassReader::open()
{
    cfp = fopen(cfname, "r");
    return true;
}

void ClassReader::getRecord()
{
    if (m_reclen == 0) return;
    size_t len = fread(buffer, sizeof(int), m_reclen, cfp);
    if (feof(cfp)) return;

    if (len != m_reclen) {
        warning("failed to read record of size %d words (%d)\n", m_reclen, len);
    }
}

int ClassReader::getInt()
{
    int data = 0;
    memcpy(&data, m_ptr, sizeof(int));
    m_ptr += sizeof(int);

    return data;
}

long int ClassReader::getLong()
{
    long int data = 0;
    memcpy(&data, m_ptr, sizeof(long int));
    m_ptr += sizeof(long int);

    return data;
}

float ClassReader::getFloat()
{
    float data = 0;
    memcpy(&data, m_ptr, sizeof(float));
    m_ptr += sizeof(float);

    return data;
}

double ClassReader::getDouble()
{
    double data = 0;
    memcpy(&data, m_ptr, sizeof(double));
    m_ptr += sizeof(double);

    return data;
}

void ClassReader::dumpRecord()
{
    int word;
    const int wpl = 8;

    if (m_reclen == 0) return;
    m_ptr = buffer;
    for (unsigned int i = 0; i < m_reclen; i++) {
        word = getInt();
        Rprintf("[%03d] %10d ", i, word);
        if ((i % wpl) == (wpl-1)) Rprintf("\n");
    }
}

void ClassReader::trim(unsigned char *input)
{
    unsigned char *dst = input, *src = input;
    unsigned char *end;

    // Skip whitespace at front...
    while (isspace(*src)) {
        ++src;
    }
    // Trim at end...
    end = src + strlen((const char *)src) - 1;
    while (end > src && isspace(*end)) {
        *end-- = 0;
    }
    // Move if needed.
    if (src != dst) {
        while ((*dst++ = *src++));
    }
}

double ClassReader::obssecond(int mjdn, double utc)
{
    const long int Jan01of1970 = 40587;
    double elapsed = 86400.0*(mjdn - Jan01of1970);
    elapsed += utc*3600.0*12.0/M_PI;

    return elapsed;
}

double ClassReader::rta(float rad)
{
    const double RADTOSEC = (3600.0*180.0/M_PI);
    return (double)rad * RADTOSEC;
}

SEXP ClassReader::headRow(int count, int scan, const char *source, const char *mol,
                          double lam, double bet,
                          double LO, double restf, double fres, double voff,
                          double time, double T, double datetime)
{
    SEXP head = PROTECT(allocVector(VECSXP, 13));
    SEXP nam = PROTECT(allocVector(STRSXP, 13)); // names attribute (column names)
    SET_STRING_ELT(nam, 0, mkChar("id"));
    SET_STRING_ELT(nam, 1, mkChar("scan"));
    SET_STRING_ELT(nam, 2, mkChar("target"));
    SET_STRING_ELT(nam, 3, mkChar("line"));
    SET_STRING_ELT(nam, 4, mkChar("RA"));
    SET_STRING_ELT(nam, 5, mkChar("Dec"));
    SET_STRING_ELT(nam, 6, mkChar("f.LO"));
    SET_STRING_ELT(nam, 7, mkChar("f0"));
    SET_STRING_ELT(nam, 8, mkChar("df"));
    SET_STRING_ELT(nam, 9, mkChar("v.LSR"));
    SET_STRING_ELT(nam,10, mkChar("dt"));
    SET_STRING_ELT(nam,11, mkChar("T.sys"));
    SET_STRING_ELT(nam,12, mkChar("observed.date"));
    namesgets(head, nam);
    UNPROTECT(1); // nam

    SEXP id = PROTECT(allocVector(INTSXP, 1));         // 0
    SEXP scanno = PROTECT(allocVector(INTSXP, 1));     // 1
    SEXP target = PROTECT(allocVector(STRSXP, 1));     // 2
    SEXP line = PROTECT(allocVector(STRSXP, 1));       // 3
    SEXP RA = PROTECT(allocVector(REALSXP, 1));        // 4
    SEXP Dec = PROTECT(allocVector(REALSXP, 1));       // 5
    SEXP fLO = PROTECT(allocVector(REALSXP, 1));       // 6
    SEXP f0 = PROTECT(allocVector(REALSXP, 1));        // 7
    SEXP df = PROTECT(allocVector(REALSXP, 1));        // 8
    SEXP vs = PROTECT(allocVector(REALSXP, 1));        // 9
    SEXP dt = PROTECT(allocVector(REALSXP, 1));        // 10
    SEXP tsys = PROTECT(allocVector(REALSXP, 1));      // 11
    // SEXP utc = PROTECT(allocVector(STRSXP, 1));        // 12
    SEXP utc = PROTECT(allocVector(REALSXP, 1));        // 12
    makePOSIXct(utc);

    INTEGER(id)[0] = count;
    INTEGER(scanno)[0] = scan;
    SET_STRING_ELT(target, 0, mkChar(source));
    SET_STRING_ELT(line, 0, mkChar(mol));
    REAL(RA)[0] = lam;
    REAL(Dec)[0] = bet;
    REAL(fLO)[0] = LO;
    REAL(f0)[0] = restf;
    REAL(df)[0] = fres;
    REAL(vs)[0] = voff;
    REAL(dt)[0] = time;
    REAL(tsys)[0] = T;
    REAL(utc)[0] = datetime;
    // SET_STRING_ELT(utc, 0, mkChar(datetime));
    SET_VECTOR_ELT(head, 0, id);
    SET_VECTOR_ELT(head, 1, scanno);
    SET_VECTOR_ELT(head, 2, target);
    SET_VECTOR_ELT(head, 3, line);
    SET_VECTOR_ELT(head, 4, RA);
    SET_VECTOR_ELT(head, 5, Dec);
    SET_VECTOR_ELT(head, 6, fLO);
    SET_VECTOR_ELT(head, 7, f0);
    SET_VECTOR_ELT(head, 8, df);
    SET_VECTOR_ELT(head, 9, vs);
    SET_VECTOR_ELT(head,10, dt);
    SET_VECTOR_ELT(head,11, tsys);
    SET_VECTOR_ELT(head,12, utc);
    UNPROTECT(13); // id, scanno, target, line, RA, Dec, fLO, f0, df, vs, dt, tsys, utc
    UNPROTECT(1);  // head
    return head;
}

SEXP ClassReader::dataVector(int nchan, float *s)
{
    m_ptr = (char *)s;
    SEXP data = PROTECT(allocVector(REALSXP, nchan));
    for (int k = 0; k < nchan; k++) REAL(data)[k] = getFloat();
    UNPROTECT(1); // data

    return data;
}

SEXP ClassReader::freqVector(int nchan, double f0, double ref, double df)
{
    SEXP freq = PROTECT(allocVector(REALSXP, nchan));
    if (df == 0.0) {
        for (int k = 0; k < nchan; k++) REAL(freq)[k] = (double)(k+1);
    } else {
        for (int k = 0; k < nchan; k++) REAL(freq)[k] = ((double)(k+1)-ref)*df + f0;
    }
    UNPROTECT(1); // freq

    return freq;
}

void ClassReader::getChar(unsigned char *dst, int len)
{
    memcpy(dst, m_ptr, len);
    dst[len] = '\0';
    trim(dst);
    m_ptr += len;
}

void ClassReader::fillHeader(char *obsblock, int code, int addr, int len)
{
    int noff = (addr-1)*4;
    m_ptr = obsblock + noff;
    if (code == -2) { /* General */
        cdesc.ut = getDouble();
        cdesc.st = getDouble();
        cdesc.az = getFloat();
        cdesc.el = getFloat();
        cdesc.tau = getFloat();
        cdesc.tsys = getFloat();
        cdesc.time = getFloat();

        if (len > 9) cdesc.xunit = getInt();
        else         cdesc.xunit = 0;
#ifdef DEBUG
        Rprintf("    -2  UT=%f %f   Az,El=%f,%f  time=%f\n", cdesc.ut, cdesc.st, cdesc.az, cdesc.el, cdesc.time);
#endif
    } else if (code == -3) { /* Position */
        if (len == 17) {
            getChar(cdesc.source, 12);
            cdesc.epoch   = getFloat();
            cdesc.lam     = getDouble();
            cdesc.bet     = getDouble();
            cdesc.lamof   = getFloat();
            cdesc.betof   = getFloat();
            cdesc.proj    = getInt();
            cdesc.sl0p    = getDouble();
            cdesc.sb0p    = getDouble();
            cdesc.sk0p    = getDouble();
        } else {
            getChar(cdesc.source, 12);
            cdesc.system  = getInt();
            cdesc.epoch   = getFloat();
            cdesc.proj    = getInt();
            cdesc.lam     = getDouble();
            cdesc.bet     = getDouble();
            cdesc.projang = getDouble();
            cdesc.lamof   = getFloat();
            cdesc.betof   = getFloat();
        }

#ifdef DEBUG
        Rprintf("    -3  '%s' %f Coord:%f,%f (%f,%f) %4d\n", cdesc.source,
               cdesc.epoch, cdesc.lam, cdesc.bet, cdesc.lamof, cdesc.betof, cdesc.proj);
#endif
    } else if (code == -4) { /* Spectroscopic */
        getChar(cdesc.line, 12);
        cdesc.restf   = getDouble();
        cdesc.nchan   = getInt();
        cdesc.rchan   = getFloat();
        cdesc.fres    = getFloat();
        cdesc.foff    = getFloat();
        cdesc.vres    = getFloat();
        cdesc.voff    = getFloat();
        cdesc.badl    = getFloat();
        cdesc.image   = getDouble();
        cdesc.vtype   = getInt();
        cdesc.doppler = getDouble();

#ifdef DEBUG
        Rprintf("    -4  '%s' %f i%f (%f %f) n=%d r=%f v=%f %f %f\n", cdesc.line,
               cdesc.restf, cdesc.image, cdesc.fres, cdesc.foff, cdesc.nchan, cdesc.rchan, cdesc.vres, cdesc.voff, cdesc.doppler);
#endif
    } else if (code == -5) { /* Baseline info section (for spectra or drifts) */
        /* Not supported, silently ignored */
    } else if (code == -6) { /* Scan numbers of initial observations */
        /* Not supported , silently ignored*/
    } else if (code == -7) { /* Default plotting limits */
        /* Not supported, silently ignored */
    } else if (code == -8) { /* Switching information */
        /* Not supported, silently ignored */
    } else if (code == -9) { /* Gauss fit results */
        /* Not supported, silently ignored */
    } else if (code == -10) { /* Continuum drifts */
        cdesc.freq     = getDouble();
        cdesc.width    = getFloat();
        cdesc.npoin    = getInt();
        cdesc.rpoin    = getFloat();
        cdesc.tref     = getFloat();
        cdesc.aref     = getFloat();
        cdesc.apos     = getFloat();
        cdesc.tres     = getFloat();
        cdesc.ares     = getFloat();
        cdesc.badc     = getFloat();
        cdesc.ctype    = getInt();
        cdesc.cimag    = getDouble();
        cdesc.colla    = getFloat();
        cdesc.colle    = getFloat();

#ifdef DEBUG
        Rprintf("   -10  f=%f %f  n=%d r=%f\n", cdesc.freq, cdesc.width, cdesc.npoin, cdesc.rpoin);
        Rprintf("   -10  t=%f %f  a=%f %f\n", cdesc.tref, cdesc.tres, cdesc.aref, cdesc.ares);
#endif
    } else if (code == -14) { /* Calibration */
        cdesc.beeff    = getFloat();
        cdesc.foeff    = getFloat();
        cdesc.gaini    = getFloat();
        cdesc.h2omm    = getFloat();
        cdesc.pamb     = getFloat();
        cdesc.tamb     = getFloat();
        cdesc.tatms    = getFloat();
        cdesc.tchop    = getFloat();
        cdesc.tcold    = getFloat();
        cdesc.taus     = getFloat();
        cdesc.taui     = getFloat();
        cdesc.tatmi    = getFloat();
        cdesc.trec     = getFloat();
        cdesc.cmode    = getInt();
        cdesc.atfac    = getFloat();
        cdesc.alti     = getFloat();
        cdesc.count[0] = getFloat();
        cdesc.count[1] = getFloat();
        cdesc.count[2] = getFloat();
        cdesc.lcalof   = getFloat();
        cdesc.bcalof   = getFloat();
        cdesc.geolong  = getDouble();
        cdesc.geolat   = getDouble();

#ifdef DEBUG
        Rprintf("   -14 %f %f %f trec=%f  (%f,%f) Site:(%f,%f, %f)\n",
               cdesc.h2omm, cdesc.tamb, cdesc.pamb, cdesc.trec,
               cdesc.lcalof, cdesc.bcalof, cdesc.geolong, cdesc.geolat, cdesc.alti);
#endif
    } else {
        warning("cannot handle CLASS section code %d yet.Sorry.", code);
    }
}

Type1Reader::Type1Reader(const char *filename) : ClassReader(filename)
{
}

Type1Reader::~Type1Reader()
{
}

void Type1Reader::getFileDescriptor()
{
    m_type = 1;

    m_reclen = 128;

    fseek(cfp, 0L, SEEK_SET);
    getRecord();
    m_ptr = buffer;

    getChar(fdesc.code, 4);

#ifdef DEBUG
    Rprintf("type=%d, reclen=%d\n", m_type, m_reclen);
#endif

    fdesc.next  = getInt();
    fdesc.lex   = getInt();
    fdesc.nex   = getInt();
    fdesc.xnext = getInt();

    int nex = fdesc.nex;
    if (nex > MAXEXT) {
        warning("number of extensions too large!");
    }
    for (int i = 0; i < nex; i++) {
        ext[i] = getInt();
#ifdef DEBUG
        Rprintf("ext[%d] = %d\n", i, ext[i]);
#endif
    }
    return;
}

void  Type1Reader::getEntry(int k)
{
    m_ptr = buffer + k*m_reclen;
    centry.xblock = getInt();
    centry.xnum = getInt();
    centry.xver = getInt();
    getChar(centry.xsourc, 12);
    getChar(centry.xline,  12);
    getChar(centry.xtel,   12);

    centry.xdobs = getInt();
    centry.xdred = getInt();
    centry.xoff1 = getFloat();
    centry.xoff2 = getFloat();
    getChar(centry.xtype, 4);
    centry.xkind = getInt();
    centry.xqual = getInt();
    centry.xscan = getInt();
    centry.xposa = getInt();
}

int Type1Reader::getDirectory()
{
    int nst = fdesc.xnext;
    long pos = (ext[0]-1)*m_reclen;

    fseek(cfp, 4*pos, SEEK_SET);

    int nrec = 2;
    int nspec = 0;
    while (nrec < nst) {
#ifdef DEBUG
        Rprintf("in %s extension 0 spectrum %d at position %ld\n", __FUNCTION__, nspec+1, ftell(cfp));
#endif
        getRecord();
        for (int k = 0; k < 4; k++) {
            getEntry(k);
            if (centry.xnum > 0 && centry.xnum < nst) {
                nspec++;
#ifdef DEBUG
                Rprintf("%5d: xblock=%6d xnum=%4d xver=%d: xsource='%s'  xline='%s' xtel='%s'\n",
                       nspec, centry.xblock, centry.xnum, centry.xver,
                       centry.xsourc, centry.xline, centry.xtel);
#endif
            } else {
                nrec = nst;
                break;
            }
        }
        nrec++;
    }
    m_nspec = nspec;
    return nspec;
}

SEXP Type1Reader::getSpectrum(int scan, bool headerOnly)
{
    long pos = (ext[0]-1)*m_reclen;
    int nrec = (scan-1)/4;
    int k = (scan-1)%4;

    fseek(cfp, 4*(pos+nrec*m_reclen), SEEK_SET);
    getRecord();
    getEntry(k);

    pos = (centry.xblock-1)*m_reclen;
    fseek(cfp, 4*pos, SEEK_SET);
    getRecord();
    m_ptr = buffer;

    getChar(csect.ident, 4);
    csect.nbl = getInt();
    csect.bytes = getInt();
    csect.adr = getInt();
    csect.nhead = getInt();
    csect.len = getInt();
    csect.ientry = getInt();
    csect.nsec = getInt();
    csect.obsnum = getInt();
    int nsec = csect.nsec;
    if (nsec > 4) nsec = 4;
    for (int i = 0; i < nsec; i++) csect.sec_cod[i] = getInt();
    for (int i = 0; i < nsec; i++) csect.sec_len[i] = getInt();
    for (int i = 0; i < nsec; i++) csect.sec_adr[i] = getInt();
#ifdef DEBUG
    Rprintf("%12s %12s %d %d %d %d %d %d %d %d %d (%d %d %d %d | %d %d %d %d | %d %d %d %d)\n",
           centry.xsourc, centry.xline, centry.xkind,
           csect.nbl, csect.bytes, csect.adr, csect.nhead, csect.len, csect.ientry, csect.nsec, csect.obsnum,
           csect.sec_cod[0], csect.sec_cod[1], csect.sec_cod[2], csect.sec_cod[3],
           csect.sec_adr[0], csect.sec_adr[1], csect.sec_adr[2], csect.sec_adr[3],
           csect.sec_len[0], csect.sec_len[1], csect.sec_len[2], csect.sec_len[3]);
#endif
    unsigned int size = csect.nbl*m_reclen*4;
    if (size > sizeof(buffer)) {
        warning("buffer too small!");
    }
    if (csect.nbl > 1) {
        size -= 4*m_reclen;
        size_t len = fread(buffer+4*m_reclen, sizeof(char), size, cfp);
        if (len != size) {
            warning("failed to read obsblock (%ld != %ld)!", len, size);
        }
    }

    for (int i = 0; i < nsec; i++) {
        fillHeader(buffer, csect.sec_cod[i], csect.sec_adr[i], csect.sec_len[i]);
    }

    bool spectrum = (centry.xkind == 0);

    double restf, rchan, LO, fres;
    if (spectrum) {
        rchan = cdesc.rchan;
        restf = cdesc.restf;
        LO = (cdesc.restf + cdesc.image)/2.0;
        fres = cdesc.fres;
    } else {
        rchan = cdesc.rpoin;
        restf = cdesc.tref;
        LO = (cdesc.freq + cdesc.cimag)/2.0;
        fres = cdesc.tres;
    }
    double lam = cdesc.lam;
    double bet = cdesc.bet;
    lam += cdesc.lamof/cos(bet);
    bet += cdesc.betof;
    double datetime = obssecond(centry.xdobs + 60549, cdesc.ut);
    // const char *datetime = obstime(centry.xdobs + 60549, cdesc.ut);
    SEXP head = PROTECT(headRow(scan, centry.xscan,
                                (const char *)centry.xsourc, (const char *)centry.xline,
                                lam*180.0/M_PI, bet*180.0/M_PI, LO, restf, fres,
                                cdesc.voff, cdesc.time, cdesc.tsys, datetime));
    if (headerOnly) {
        UNPROTECT(1); // head
        return head;
    }

    int ndata = 0;
    if (spectrum) ndata = cdesc.nchan;
    else          ndata = cdesc.npoin;

    if (ndata > (int)MAXCHANNELS) {
        Rprintf("maximum number of channels exceeded: %d %d\n", ndata, MAXCHANNELS);
    }

    float *s = (float *)(buffer+4*(csect.nhead-1));
    SEXP data = PROTECT(dataVector(ndata, s));
    SEXP freq = PROTECT(freqVector(ndata, restf, rchan, fres));

    SEXP S = PROTECT(makeSpectrum(head, freq, data));
    UNPROTECT(4); // head, freq, data, S

    return S;
}

Type2Reader::Type2Reader(const char *filename) : ClassReader(filename)
{
}

Type2Reader::~Type2Reader()
{
}

void Type2Reader::getFileDescriptor()
{
    size_t len = 0;
    m_type = 2;

    fseek(cfp, 0L, SEEK_SET);
    len = fread((char *)&fdesc, sizeof(char), 4, cfp);
    len = fread((char *)&fdesc.reclen, sizeof(int), 1, cfp);
    if (len != 1) {
        warning("failed to read record length\n");
    }
    m_reclen = fdesc.reclen;

#ifdef DEBUG
    Rprintf("type=%d, reclen=%d\n", m_type, m_reclen);
#endif

    fseek(cfp, 0L, SEEK_SET);
    getRecord();
    m_ptr = buffer+4*sizeof(char)+sizeof(int);

    fdesc.kind  = getInt();
    fdesc.vind  = getInt();
    fdesc.lind  = getInt();
    fdesc.flags = getInt();

    fdesc.xnext    = getLong();
    fdesc.nextrec  = getLong();
    fdesc.nextword = getInt();

    fdesc.lex1 = getInt();
    fdesc.nex  = getInt();
    fdesc.gex  = getInt();

    if (fdesc.kind != 1) {
        warning("not a file written by CLASS!");
        return;
    }
    if ((fdesc.gex != 10) && (fdesc.gex != 20)) {
        warning("problem with extension growth!");
        return;
    }

    int nex = fdesc.nex;
    if (nex > MAXEXT) {
        warning("number of extensions too large!");
    }
    unsigned int size = m_reclen;
    for (int i = 0; i < nex; i++) {
        ext[i] = getLong();
        if (fdesc.gex == 20) {
            for (int j = 0; j < i; j++) size *= 2;
        }
#ifdef DEBUG
        Rprintf("ext[%d] = %ld %d\n", i, ext[i], size);
#endif
    }
    return;
}

void Type2Reader::getEntry(int k)
{
    m_ptr = buffer + k*4*fdesc.lind;
    centry.xblock = getLong();
    centry.xword = getInt();
    centry.xnum = getLong();
    centry.xver = getInt();
    getChar(centry.xsourc, 12);
    getChar(centry.xline,  12);
    getChar(centry.xtel,   12);

    centry.xdobs = getInt();
    centry.xdred = getInt();
    centry.xoff1 = getFloat();
    centry.xoff2 = getFloat();
    getChar(centry.xtype, 4);
    centry.xkind = getInt();
    centry.xqual = getInt();
    centry.xposa = getInt();
    centry.xscan = getLong();
    centry.xsubs = getInt();
}

int Type2Reader::getDirectory()
{
    int growth = 1;

    int nspec = 0;
    for (int iext = 0; iext < fdesc.nex; iext++) {
        int nst = fdesc.lex1*growth;
        unsigned int isize = nst*fdesc.lind;
        long pos = (ext[iext]-1)*1024;

        fseek(cfp, 4*pos, SEEK_SET);
#ifdef DEBUG
        Rprintf("in %s extension %d spectrum %d at position %ld\n", __FUNCTION__, iext, nspec+1, ftell(cfp));
#endif
        size_t len = fread(buffer, sizeof(int), isize, cfp);
        if (len != isize) {
            warning("failed to read directory information\n");
            return 0;
        }
        m_ptr = buffer;

        for (int k = 0; k < nst; k++) {
            getEntry(k);
            if (centry.xnum >= 1) {
                nspec++;
#ifdef DEBUG
                Rprintf("%5d: xblock=%3ld xword=%4d xnum=%3ld xver=%d: xsource='%s'  xline='%s' xtel='%s'\n",
                       nspec, centry.xblock, centry.xword, centry.xnum, centry.xver,
                       centry.xsourc, centry.xline, centry.xtel);
#endif
            }
        }
        if (fdesc.gex == 20) growth *= 2;
    }
    m_nspec = nspec;
    return nspec;
}

SEXP Type2Reader::getSpectrum(int scan, bool headerOnly)
{
    size_t len = 0;
    int iext = 0;
    unsigned int isize = 0;
    int jent = 0;

    int growth = 1;
    int nspec = 0;
    for (int i = 0; i < fdesc.nex; i++) {
        int nst = fdesc.lex1*growth;
        isize = nst*fdesc.lind;
        for (int k = 0; k < nst; k++) {
            nspec++;
            if (nspec == scan) {
                jent = k;
                iext = i;
                break;
            }
        }
        if (nspec == scan) break;
        if (fdesc.gex == 20) growth *= 2;
    }

    long pos = (ext[iext]-1)*m_reclen;
#ifdef DEBUG
    Rprintf("in %s spectrum %d at position %ld in extension %d entry %d at word %d\n",
            __FUNCTION__, scan, pos, iext, jent, startword);
#endif
    fseek(cfp, 4*pos, SEEK_SET);
    len = fread(buffer, sizeof(int), isize, cfp);
    if (len != isize) {
        warning("failed to read entry descriptor\n");
        return R_NilValue;
    }

    getEntry(jent);
    pos = (centry.xblock-1)*m_reclen+centry.xword-1;

    fseek(cfp, 4*pos, SEEK_SET);
#ifdef DEBUG
    Rprintf("in %s reading spectrum %d at block %ld, word %d position %ld (%ld)\n",
            __FUNCTION__, scan, centry.xblock, centry.xword, ftell(cfp), 4*pos);
#endif
    getRecord();
    m_ptr = buffer;

    getChar(csect.ident, 4);
    csect.version = getInt();
    csect.nsec = getInt();
    csect.nword = getLong();
    csect.adata = getLong();
    csect.ldata = getLong();
    csect.xnum = getLong();
    int nsec = csect.nsec;
    if (nsec > 10) nsec = 10;
    for (int i = 0; i < nsec; i++) csect.sec_cod[i] = getInt();
    for (int i = 0; i < nsec; i++) csect.sec_len[i] = getLong();
    for (int i = 0; i < nsec; i++) csect.sec_adr[i] = getLong();
#ifdef DEBUG
    Rprintf("%7d %12s %12s %d %d (%d %d %d %d | %d %d %d %d | %d %d %d %d)\n",
           scan,
           centry.xsourc, centry.xline, centry.xkind,
           csect.nsec,
           csect.sec_cod[0], csect.sec_cod[1], csect.sec_cod[2], csect.sec_cod[3],
           csect.sec_adr[0], csect.sec_adr[1], csect.sec_adr[2], csect.sec_adr[3],
           csect.sec_len[0], csect.sec_len[1], csect.sec_len[2], csect.sec_len[3]);
#endif

    for (int i = 0; i < nsec; i++) {
        unsigned int isize = csect.sec_len[i];
        unsigned int size = 4*isize;
        long secpos = 4*(pos + csect.sec_adr[i]-1);
#ifdef DEBUG
        Rprintf("in %s reading section %d at position %ld (%ld,%ld)\n",
                __FUNCTION__, i, ftell(cfp), pos, csect.sec_adr[i]);
#endif
        int ierr = fseek(cfp, secpos, SEEK_SET);
        if (ierr == -1) {
            warning("failure to position section %d (%ld)\n", i, secpos);
            break;
        }
        len = fread(buffer, sizeof(char), size, cfp);
        if (len != size) {
            warning("failed to read section %d (%d != %d)\n", i, size, len);
            break;
        }
        m_ptr = buffer;
#ifdef DEBUG
        Rprintf("Sect[%d]=%3d at adr %d of length %d.\n", i, csect.sec_cod[i], secpos, size);
#endif
        fillHeader(buffer, csect.sec_cod[i], 1, csect.sec_len[i]);
    }

    bool spectrum = (centry.xkind == 0);

    double restf, rchan, LO, fres;
    if (spectrum) {
        rchan = cdesc.rchan;
        restf = cdesc.restf;
        LO = (cdesc.restf + cdesc.image)/2.0;
        fres = cdesc.fres;
    } else {
        rchan = cdesc.rpoin;
        restf = cdesc.tref;
        LO = (cdesc.freq + cdesc.cimag)/2.0;
        fres = cdesc.tres;
    }
    double lam = cdesc.lam;
    double bet = cdesc.bet;
    lam += cdesc.lamof/cos(bet);
    bet += cdesc.betof;
    double datetime = obssecond(centry.xdobs + 60549, cdesc.ut);
    // const char *datetime = obstime(centry.xdobs + 60549, cdesc.ut);
    SEXP head = PROTECT(headRow(scan, centry.xscan,
                                (const char *)centry.xsourc, (const char *)centry.xline,
                                lam*180.0/M_PI, bet*180.0/M_PI, LO, restf, fres,
                                cdesc.voff, cdesc.time, cdesc.tsys, datetime));
    if (headerOnly) {
        UNPROTECT(1); // head
        return head;
    }

    isize = csect.ldata;
    char *datatbl = buffer;
    if (4*isize > sizeof(buffer)) {
        warning("buffer too small!");
    }
    long datapos = 4*(pos + csect.adata-1);
    fseek(cfp, datapos, SEEK_SET);
#ifdef DEBUG
    Rprintf("in %s reading data for spectrum %d at position %ld (%ld,%ld)\n",
            __FUNCTION__, scan, ftell(cfp), pos, csect.adata);
#endif
    len = fread(datatbl, sizeof(int), isize, cfp);
    if (len != isize) {
        warning("failed to read data block\n");
    }

    int ndata = 0;
    if (spectrum) ndata = cdesc.nchan;
    else          ndata = cdesc.npoin;

    if (ndata > (int)MAXCHANNELS) {
        Rprintf("maximum number of channels exceeded: %d %d\n", ndata, MAXCHANNELS);
    }

    float *s = (float *)(datatbl);
    SEXP data = PROTECT(dataVector(ndata, s));
    SEXP freq = PROTECT(freqVector(ndata, restf, rchan, fres));

    SEXP S = PROTECT(makeSpectrum(head, freq, data));
    UNPROTECT(4); // head, freq, data, S

    return S;
}
