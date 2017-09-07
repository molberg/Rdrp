#ifndef CLASSREADER_H
#define CLASSREADER_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

int fileType(const char *filename);
int qcmp(const void *x, const void *y);
SEXP makeFactor(SEXP str);
SEXP emptyFrame(int nspec);
SEXP makeSpectrum(SEXP head, SEXP freq, SEXP data);
SEXP getClassHeader(SEXP filename);
SEXP readClass(SEXP filename, SEXP header = R_NilValue);

struct ClassDescriptor {
    int xbloc;
    int xnum;
    int xver;
    unsigned char xsourc[13];
    char xline[13];
    char xtel[13];
    int xdobs;
    int xdred;
    float xoff1;
    float xoff2;
    char xtype[5];
    int xkind;
    int xqual;
    int xscan;
    int xposa;
    /* Version 2 additions */
    long long int xlbloc, xlnum, xlscan;
    int xword, xsubs;
    /* End */
    char xfront[9];
    char xback[9];
    char xproc[9];
    char xproj[9];
    char unused[25];

    double ut, st;
    float az, el;
    float tau, tsys, time;
    int xunit;

    unsigned char source[13];
    float epoch;
    double lam, bet, projang;
    float lamof, betof;
    int system, proj;
    double sl0p, sb0p, sk0p;

    unsigned char line[13];
    double restf;
    int nchan;
    float rchan, fres, foff, vres, voff, badl;
    double image;
    int vtype;
    double doppler;

    float beeff, foeff, gaini, h2omm, pamb, tamb, tatms, tchop, tcold, taus, taui, tatmi, trec;
    int cmode;
    float atfac, alti, count[3], lcalof, bcalof;
    double geolong, geolat;

    double freq;
    float width;
    int npoin;
    float rpoin, tref, aref, apos, tres, ares, badc;
    int ctype;
    double cimag;
    float colla, colle;

    int ndata;
    float *data;
};

struct FileDescriptor1 {
    unsigned char code[4];
    int next;
    int lex;
    int nex;
    int lind;
    int xnext;
};

struct FileDescriptor2 {
    unsigned char code[4];
    int reclen;
    int kind;
    int vind;
    int lind;
    int flags;
    long int xnext;
    long int nextrec;
    int nextword;
    int lex1;
    int nex;
    int gex;
};

struct Type1Entry {
    int xblock;
    int xnum;
    int xver;
    unsigned char xsourc[13];
    unsigned char xline[13];
    unsigned char xtel[13];
    unsigned char xtype[5];
    int xdobs;
    int xdred;
    float xoff1;
    float xoff2;
    int xkind;
    int xqual;
    int xposa;
    int xscan;
    unsigned char xfront[9];
    unsigned char xback[9];
    unsigned char xproc[9];
    unsigned char xproj[9];
};

struct Type2Entry {
    long int xblock;
    int xword;
    long int xnum;
    int xver;
    unsigned char xsourc[13];
    unsigned char xline[13];
    unsigned char xtel[13];
    unsigned char xtype[5];
    int xdobs;
    int xdred;
    float xoff1;
    float xoff2;
    int xkind;
    int xqual;
    int xposa;
    long int xscan;
    int xsubs;
};

struct ClassSection1 {
    unsigned char ident[5];
    int nbl;
    int bytes;
    int adr;
    int nhead;
    int len;
    int ientry;
    int nsec;
    int obsnum;
    int sec_cod[4];
    int sec_adr[4];
    int sec_len[4];
};

struct ClassSection2 {
    unsigned char ident[5];
    int version;
    int nsec;
    long int nword;
    long int adata;
    long int ldata;
    long int xnum;
    int sec_cod[10];
    long int sec_len[10];
    long int sec_adr[10];
};

class ClassReader {
 public:
    ClassReader(const char *);
    virtual ~ClassReader();

    bool open();
    virtual void getFileDescriptor() = 0;
    virtual int  getDirectory() = 0;
    virtual SEXP getSpectrum(int scan, bool headerOnly = false) = 0;
    void dumpRecord();
    
 protected:
    void getRecord();
    int getInt();
    long int getLong();
    void getChar(unsigned char *dst, int len);
    float getFloat();
    double getDouble();
    void fillHeader(char *obsblock, int code, int addr, int len);
    void fillData(char *obsblock, int nhead, int ndata);
    void trim(unsigned char *ptr);
    const char *obstime(int mjdn, double utc);
    double rta(float rad);
    SEXP headRow(int count, int scan,
                 const char *source, const char *mol,
                 double lam, double bet,
                 double LO, double restf, double fres,
                 double voff, double time, double T,
                 const char *datetime);
    SEXP dataVector(int nchan, float *data);
    SEXP freqVector(int nchan, double f0, double ref, double df);

    struct ClassDescriptor cdesc;
    float *m_data;
    int m_type;
    char cfname[256];
    FILE *cfp;
    char *record;
    char *dir;
    char *m_ptr;
    int m_reclen;
    int m_nspec;
};

class Type1Reader : public ClassReader {

 public:
    Type1Reader(const char *);
    ~Type1Reader();

    void getFileDescriptor();
    int  getDirectory();
    SEXP getSpectrum(int scan, bool headerOnly = false);

 private:
    FileDescriptor1 fdesc;
    Type1Entry centry;
    ClassSection1 csect;
    int *ext;
};

class Type2Reader : public ClassReader {

 public:
    Type2Reader(const char *);
    ~Type2Reader();

    void getFileDescriptor();
    int  getDirectory();
    SEXP getSpectrum(int scan, bool headerOnly = false);

 private:
    FileDescriptor2 fdesc;
    Type2Entry centry;
    ClassSection2 csect;
    long int *ext;
};

#endif
