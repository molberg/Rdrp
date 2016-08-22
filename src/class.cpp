#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "class.h"

#define RADTOSEC  (3600.0*180.0/M_PI)

#undef DEBUG
#undef ALLOC

const char *obstime(int mjdn, double utc)
{
    static char datetime[32];
    long int jdn, j, g, dg, c, dc, b, db, a, da, y, m , d;
    int year, month, day, hour, min, sec;
    double ip;

    jdn =  mjdn + 2400001;

    j = jdn + 32044;
    g = j / 146097;          dg = j % 146097;
    c = dg / 36524;          dc = dg - c * 36524;
    b = dc / 1461;           db = dc % 1461;
    a = ((db/365 + 1)*3)/4 ; da = db - a*365;
    y = g*400 + c*100 + b*4 + a;
    m = (da*5 + 308)/153 - 2;
    d = da - ((m+4)*153)/5 + 122;
    year = y - 4800 + (m+2)/12;
    month = ((m+2) % 12) + 1;
    day = d + 1;
    utc *= 12.0/M_PI;
    utc = modf(utc, &ip);
    hour  = (int)ip;
    utc = modf(60.0*utc, &ip);
    min = (int)ip;
    utc = modf(60.0*utc, &ip);
    sec = (int)ip;

    sprintf(datetime, "%04d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, min, sec);
    return datetime;
}

double rta(float rad)
{
    return (double)rad * RADTOSEC;
}

int nst = 0, current = 0;
CLASS *st = NULL;

void trim(char *input)
{
    char *dst = input, *src = input;
    char *end;

    // Skip whitespace at front...
    while (isspace((unsigned char)*src)) {
        ++src;
    }
    // Trim at end...
    end = src + strlen(src) - 1;
    while (end > src && isspace((unsigned char)*end)) {
        *end-- = 0;
    }
    // Move if needed.
    if (src != dst) {
        while ((*dst++ = *src++));
    }
}

void fill_class_v1_obs(CLASS_SECTION *o, char *edesc)
{
    int k, nmax, n;

    if (!o) return;
    if (!edesc) return;

    n = 0;
    memcpy(o->ident,     &edesc[n], 4); n += 4;
    memcpy(&(o->nbl),    &edesc[n], 4); n += 4;
    memcpy(&(o->bytes),  &edesc[n], 4); n += 4;
    memcpy(&(o->adr),    &edesc[n], 4); n += 4;
    memcpy(&(o->nhead),  &edesc[n], 4); n += 4;
    memcpy(&(o->len),    &edesc[n], 4); n += 4;
    memcpy(&(o->ientry), &edesc[n], 4); n += 4;
    memcpy(&(o->nsec),   &edesc[n], 4); n += 4;
    memcpy(&(o->obsnum), &edesc[n], 4); n += 4;

    nmax = o->nsec;
    if (nmax >= 4) nmax=4;
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_cod[k]), &edesc[n], 4); n += 4;
    }
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_len[k]), &edesc[n], 4); n += 4;
    }
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_adr[k]), &edesc[n], 4); n += 4;
    }
}

void fill_class_v1_table(CLASS *o, char *i)
{
    int n;

    if (!o) return;
    if (!i) return;

    n = 0;
    memcpy(&(o->xbloc), &i[n],  4); n +=  4;
    memcpy(&(o->xnum),  &i[n],  4); n +=  4;
    memcpy(&(o->xver),  &i[n],  4); n +=  4;
    strncpy(o->xsourc,  &i[n], 12); n += 12; o->xsourc[12] = '\0';
    strncpy(o->xline,   &i[n], 12); n += 12; o->xline[12]  = '\0';
    strncpy(o->xtel,    &i[n], 12); n += 12; o->xtel[12]   = '\0';
    memcpy(&(o->xdobs), &i[n],  4); n +=  4;
    memcpy(&(o->xdred), &i[n],  4); n +=  4;
    memcpy(&(o->xoff1), &i[n],  4); n +=  4;
    memcpy(&(o->xoff2), &i[n],  4); n +=  4;
    strncpy(o->xtype,   &i[n],  4); n +=  4; o->xtype[4]   = '\0';
    memcpy(&(o->xkind), &i[n],  4); n +=  4;
    memcpy(&(o->xqual), &i[n],  4); n +=  4;
    memcpy(&(o->xscan), &i[n],  4); n +=  4;
    memcpy(&(o->xposa), &i[n],  4); n +=  4;
    strncpy(o->xfront,  &i[n],  8); n +=  8; o->xfront[8]  = '\0';
    strncpy(o->xback,   &i[n],  8); n +=  8; o->xback[8]   = '\0';
    strncpy(o->xproc,   &i[n],  8); n +=  8; o->xproc[8]   = '\0';
    strncpy(o->xproj,   &i[n],  8); n +=  8; o->xproj[8]   = '\0';
    strncpy(o->unused,  &i[n], 24); n += 24; o->unused[24] = '\0';
}

void fill_class_v2_obs(CLASS_SECTION_v2 *o, char *edesc)
{
    int k, nmax, n;

    if (!o) return;
    if (!edesc) return;

    n = 0;
    memcpy(o->ident,      &edesc[n], 4); n += 4;
    memcpy(&(o->version), &edesc[n], 4); n += 4;
    memcpy(&(o->nsec),    &edesc[n], 4); n += 4;
    memcpy(&(o->nword),   &edesc[n], 8); n += 8;
    memcpy(&(o->adata),   &edesc[n], 8); n += 8;
    memcpy(&(o->ldata),   &edesc[n], 8); n += 8;
    memcpy(&(o->xnum),    &edesc[n], 8); n += 8;

    nmax = o->nsec;
    if (nmax >= MAX_CLASS_SECT) nmax = MAX_CLASS_SECT;
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_cod[k]), &edesc[n], 4); n += 4;
    }
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_len[k]), &edesc[n], 8); n += 8;
    }
    for (k = 0; k < nmax; k++) {
        memcpy(&(o->sec_adr[k]), &edesc[n], 8); n += 8;
    }
#ifdef DEBUG
    Rprintf("n=%d\n", n);
    Rprintf("  ident   '%s' Identifier of Entry Description.\n", o->ident);
    Rprintf("  version '%d' Observations version.\n", o->version);
    Rprintf("  nsec    '%d' Number of sections.\n", o->nsec);
    Rprintf("  nword   '%lli' Length of this entry [words].\n", o->nword);
    Rprintf("  adata   '%lli' Data address [word].\n", o->adata);
    Rprintf("  ldata   '%lli' Data length [words].\n", o->ldata);
    Rprintf("  xnum    '%lli' Entry number.\n", o->xnum);
    for (k = 0; k < nmax; k++) {
        Rprintf("  %3d len=%lld words at addr %lld\n", o->sec_cod[k], o->sec_len[k], o->sec_adr[k]);
    }
#endif
}

void fill_class_v2_table(CLASS *o, char *i)
{
    int n;

    if (!o) return;
    if (!i) return;

    n = 0;
    memcpy(&(o->xlbloc), &i[n], 8); n += 8;
    memcpy(&(o->xword), &i[n],  4); n += 4;
    memcpy(&(o->xlnum), &i[n],  8); n += 8;
    memcpy(&(o->xver), &i[n],   4); n += 4;
    strncpy(o->xsourc,  &i[n], 12); n += 12; o->xsourc[12] = '\0';
    strncpy(o->xline,   &i[n], 12); n += 12; o->xline[12]  = '\0';
    strncpy(o->xtel,    &i[n], 12); n += 12; o->xtel[12]   = '\0';
    memcpy(&(o->xdobs), &i[n],  4); n += 4;
    memcpy(&(o->xdred), &i[n],  4); n += 4;
    memcpy(&(o->xoff1), &i[n],  4); n += 4;
    memcpy(&(o->xoff2), &i[n],  4); n += 4;
    strncpy(o->xtype,   &i[n],  4); n +=  4; o->xtype[4]   = '\0';
    memcpy(&(o->xkind), &i[n],  4); n += 4;
    memcpy(&(o->xqual), &i[n],  4); n += 4;
    memcpy(&(o->xposa), &i[n],  4); n += 4;
    memcpy(&(o->xlscan), &i[n], 4); n += 8;
    memcpy(&(o->xsubs), &i[n],  4); n += 4;

    o->xbloc = (int)o->xlbloc;
    o->xnum = (int)o->xlnum;
    o->xscan = (int)o->xlscan;
}

void list_xdata(int nbl, CLASS *c)
{
    Rprintf("%5d 1 %7d %4d %d %12s %12s %12s %d (%+7.1f\",%+7.1f\") %d\n",
            nbl+1,
            c->xbloc, c->xnum, c->xver, c->xsourc, c->xline, c->xtel,
            c->xdobs, rta(c->xoff1), rta(c->xoff2), c->xscan);
}

int get_class_v1_file_listing(FILE *fp, CLASS_INFO *info)
{
    int n, k, len, done, nscan;
    static char block[CLASS_BLK_SIZE];
    static CLASS obs;

    if (!fp) return -1;

    nst = info->first;
    st = (CLASS *)malloc(nst*sizeof(CLASS));
    if (!st) {
        nst = 0;
        warning("Out of memory: Cannot allocate CLASS file structure.");
        fclose(fp);
        return -1;
    }
#ifdef ALLOC
    Rprintf("%p -- allocated %d CLASS structures\n", st, nst);
#endif
    obs.data = NULL;

    /* move one blk forward */
    fseek(fp, CLASS_BLK_SIZE*sizeof(char), SEEK_CUR);
    n = 2;
    done = 0;
    nscan = 0;
    do {
        len = fread(block, sizeof(char), CLASS_BLK_SIZE, fp);
        if (len != CLASS_BLK_SIZE) {
            warning("Cannot read %d byte block from CLASS file.", CLASS_BLK_SIZE);
            fclose(fp);
            return -2;
        }

        for (k = 0; k < CLASS_BLK_SIZE; k += CLASS_BLK_SIZE/4) {
            fill_class_v1_table(&obs, block+k);
            if (obs.xnum > 0 && obs.xnum < nst) {
                st[obs.xnum] = obs;
                list_xdata(n, &obs);
                nscan++;
            } else {
                done = 1;
                break;
            }
        }
        if (!done) n++;
    } while (!done && (n < nst));

    return nscan;
}

int get_class_v2_file_listing(FILE *fp, CLASS_INFO *info, int extno, int scan_only)
{
    int isize, len, j, nread, pos, extgrowth = 1;
    char *index = NULL, *iptr;
    CLASS obs;

    if (!fp) return -1;

    for (j = 0; j < extno; j++) extgrowth *= 2;
    pos = (info->aex[extno] - 1)*1024;
    isize = 4 * info->lex1 * extgrowth * info->lind;
    nst = info->first * extgrowth;

#ifdef DEBUG
    Rprintf("*** extno=%d   nst=%d  pos=%d   isize=%d\n", extno, nst, pos, isize);
#endif

    st = (CLASS *)malloc(nst*sizeof(CLASS));
    if (!st) {
        nst = 0;
        warning("Out of memory: Cannot allocate CLASS file structure.");
        fclose(fp);
        return -1;
    }
#ifdef ALLOC
    Rprintf("%p -- allocated %d CLASS structures\n", st, nst);
#endif

    /* lex1 is max no of entries in the first extension block */
    /* nex is no of extension blocks, starting at record aex[] */
    /* next available entry is xnext, so we have valid entries between [1,xnext-1] */

    /* size=4*lex1*lind should hold maximum possible index */
    index = (char *)malloc(isize * sizeof(char));
    if (!index) {
        warning("Cannot allocate CLASS index of length %d.", isize);
        fclose(fp);
        return -1;
    }
#ifdef ALLOC
    Rprintf("%p -- allocated %d characters for index\n", index, isize);
#endif

    fseek(fp, 4*pos, SEEK_SET);
    len = fread(index, sizeof(char), isize, fp);
    if (len != isize) {
        warning("Read only %d bytes of %d from CLASS file.", len, isize);
        free(index);
#ifdef ALLOC
        Rprintf("%p -- free\n", index);
#endif
        fclose(fp);
        return -2;
    }

    nread = 0;
    iptr = index;
    for (j = 0; j < nst; j++) { /* loop over all existing indices */
        fill_class_v2_table(&obs, iptr);
        obs.data = NULL;
        st[j] = obs;
        if (obs.xnum >= 1) {
            if (!scan_only) list_xdata(j, &obs);
            nread++;
        }
        /* move pointer to next entry */
        iptr += 4*info->lind;
    }
    free(index);
#ifdef ALLOC
    Rprintf("%p -- free\n", index);
#endif
    /* return the no of scans */
    return nread;
}

int check_block(CLASS *c, int nst, int block_no)
{
    int n;

    for (n = 1; n < nst; n++) {
        if (c[n].xbloc == block_no) return n;
    }

    return 0;
}

void fill_class_header(int code, int first, int len, char *s, int size, CLASS *c)
{
    int n;

    if (!s) return;
    if (!c) return;

    n = (first-1)*4;
    if (code == -2) { /* General */
        memcpy(&(c->g.ut),   &s[n], 8); n += 8;
        memcpy(&(c->g.st),   &s[n], 8); n += 8;
        memcpy(&(c->g.az),   &s[n], 4); n += 4;
        memcpy(&(c->g.el),   &s[n], 4); n += 4;
        memcpy(&(c->g.tau),  &s[n], 4); n += 4;
        memcpy(&(c->g.tsys), &s[n], 4); n += 4;
        memcpy(&(c->g.time), &s[n], 4); n += 4;

        if (len > 9) {
            memcpy(&(c->g.xunit), &s[n], 4); n += 4;
        } else {
            c->g.xunit = 0;
        }
#ifdef DEBUG
        Rprintf("    -2  UT=%f %f   Az,El=%f,%f  time=%f\n", c->g.ut, c->g.st, c->g.az, c->g.el, c->g.time);
#endif
    } else if (code == -3) { /* Position */
        strncpy(c->p.source,  &s[n], 12); n += 12; c->p.source[12] = '\0';
        memcpy(&(c->p.epoch), &s[n],  4); n +=  4;
        memcpy(&(c->p.lam),   &s[n],  8); n +=  8;
        memcpy(&(c->p.bet),   &s[n],  8); n +=  8;
        memcpy(&(c->p.lamof), &s[n],  4); n +=  4;
        memcpy(&(c->p.betof), &s[n],  4); n +=  4;
        memcpy(&(c->p.proj),  &s[n],  4); n +=  4;
        memcpy(&(c->p.sl0p),  &s[n],  8); n +=  8;
        memcpy(&(c->p.sb0p),  &s[n],  8); n +=  8;
        memcpy(&(c->p.sk0p),  &s[n],  8); n +=  8;

#ifdef DEBUG
        Rprintf("    -3  '%s' %f Coord:%f,%f (%f,%f) %4d\n", c->p.source,
                c->p.epoch, c->p.lam, c->p.bet, c->p.lamof, c->p.betof, c->p.proj);
#endif
    } else if (code == -4) { /* Spectroscopic */
        strncpy(c->s.line,      &s[n], 12); n += 12; c->s.line[12] = '\0';
        memcpy(&(c->s.restf),   &s[n],  8); n += 8;
        memcpy(&(c->s.nchan),   &s[n],  4); n += 4;
        memcpy(&(c->s.rchan),   &s[n],  4); n += 4;
        memcpy(&(c->s.fres),    &s[n],  4); n += 4;
        memcpy(&(c->s.foff),    &s[n],  4); n += 4;
        memcpy(&(c->s.vres),    &s[n],  4); n += 4;
        memcpy(&(c->s.voff),    &s[n],  4); n += 4;
        memcpy(&(c->s.bad),     &s[n],  4); n += 4;
        memcpy(&(c->s.image),   &s[n],  8); n += 8;
        memcpy(&(c->s.vtype),   &s[n],  4); n += 4;
        memcpy(&(c->s.doppler), &s[n],  8); n += 8;

#ifdef DEBUG
        Rprintf("    -4  '%s' %f i%f (%f %f) n=%d r=%f v=%f %f %f\n", c->s.line,
                c->s.restf, c->s.image, c->s.fres, c->s.foff, c->s.nchan, c->s.rchan, c->s.vres, c->s.voff, c->s.doppler);
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
        memcpy(&(c->u.freq),  &s[n], 8); n += 8;
        memcpy(&(c->u.width), &s[n], 4); n += 4;
        memcpy(&(c->u.npoin), &s[n], 4); n += 4;
        memcpy(&(c->u.rpoin), &s[n], 4); n += 4;
        memcpy(&(c->u.tref),  &s[n], 4); n += 4;
        memcpy(&(c->u.aref),  &s[n], 4); n += 4;
        memcpy(&(c->u.apos),  &s[n], 4); n += 4;
        memcpy(&(c->u.tres),  &s[n], 4); n += 4;
        memcpy(&(c->u.ares),  &s[n], 4); n += 4;
        memcpy(&(c->u.bad),   &s[n], 4); n += 4;
        memcpy(&(c->u.ctype), &s[n], 4); n += 4;
        memcpy(&(c->u.cimag), &s[n], 8); n += 8;
        memcpy(&(c->u.colla), &s[n], 4); n += 4;
        memcpy(&(c->u.colle), &s[n], 4); n += 4;

#ifdef DEBUG
        Rprintf("   -10  f=%f %f  n=%d r=%f\n", c->u.freq, c->u.width, c->u.npoin, c->u.rpoin);
        Rprintf("   -10  t=%f %f  a=%f %f\n", c->u.tref, c->u.tres, c->u.aref, c->u.ares);
#endif
    } else if (code == -14) { /* Calibration */
        memcpy(&(c->c.beeff),    &s[n], 4); n += 4;
        memcpy(&(c->c.foeff),    &s[n], 4); n += 4;
        memcpy(&(c->c.gaini),    &s[n], 4); n += 4;
        memcpy(&(c->c.h2omm),    &s[n], 4); n += 4;
        memcpy(&(c->c.pamb),     &s[n], 4); n += 4;
        memcpy(&(c->c.tamb),     &s[n], 4); n += 4;
        memcpy(&(c->c.tatms),    &s[n], 4); n += 4;
        memcpy(&(c->c.tchop),    &s[n], 4); n += 4;
        memcpy(&(c->c.tcold),    &s[n], 4); n += 4;
        memcpy(&(c->c.taus),     &s[n], 4); n += 4;
        memcpy(&(c->c.taui),     &s[n], 4); n += 4;
        memcpy(&(c->c.tatmi),    &s[n], 4); n += 4;
        memcpy(&(c->c.trec),     &s[n], 4); n += 4;
        memcpy(&(c->c.cmode),    &s[n], 4); n += 4;
        memcpy(&(c->c.atfac),    &s[n], 4); n += 4;
        memcpy(&(c->c.alti),     &s[n], 4); n += 4;
        memcpy(&(c->c.count[0]), &s[n], 4); n += 4;
        memcpy(&(c->c.count[1]), &s[n], 4); n += 4;
        memcpy(&(c->c.count[2]), &s[n], 4); n += 4;
        memcpy(&(c->c.lcalof),   &s[n], 4); n += 4;
        memcpy(&(c->c.bcalof),   &s[n], 4); n += 4;
        memcpy(&(c->c.geolong),  &s[n], 8); n += 8;
        memcpy(&(c->c.geolat),   &s[n], 8); n += 8;

#ifdef DEBUG
        Rprintf("   -14 %d %f %f %f trec=%f  (%f,%f) Site:(%f,%f, %f)\n", n/4, c->c.h2omm, c->c.tamb, c->c.pamb, c->c.trec,
                c->c.lcalof, c->c.bcalof, c->c.geolong, c->c.geolat, c->c.alti);
#endif
    } else {
        warning("Cannot handle CLASS section code %d yet.Sorry.", code);
    }
}

int fill_class_v1_data(CLASS_SECTION *cs, char *s, int size, CLASS *c)
{
    int i, n, ndata;
    char *p;

    n = 4*(cs->nhead - 1);

    if (c->xkind == 1) ndata = c->u.npoin;
    else               ndata = c->s.nchan;

    c->ndata = ndata;
    c->data = (float *)malloc(ndata * sizeof(float));
    if (!c->data) return -1;
#ifdef ALLOC
    Rprintf("%p -- allocated %d floats for data\n", c->data, ndata);
#endif
    p = &s[n];
    for (i = 0; i < ndata; i++) {
        memcpy(&(c->data[i]), p, 4);
        p += 4;
    }

#ifdef DEBUG
    if (ndata > 30)
        Rprintf("    -D %d %d  %f %f %f %f ... %f %f\n", ndata, n,
                c->data[0], c->data[10], c->data[20], c->data[30],
                c->data[ndata-2], c->data[ndata-1]);
#endif
    return ndata;
}

int fill_class_v2_data(CLASS_SECTION_v2 *cs, char *s, int size, CLASS *c)
{
    int i, ndata;
    char *p;

    if (c->xkind == 1) ndata = c->u.npoin;
    else               ndata = c->s.nchan;

    c->ndata = ndata;
    c->data = (float *)malloc(ndata * sizeof(float));
    if (!c->data) return -1;
#ifdef ALLOC
    Rprintf("%p -- allocated %d floats for data\n", c->data, ndata);
#endif
    p = &s[0];
    for (i = 0; i < ndata; i++) {
        memcpy(&(c->data[i]), p, 4);
        p += 4;
    }

#ifdef DEBUG
    if (ndata > 30)
        Rprintf("    -D2 %d  %f %f %f %f ... %f %f\n", ndata,
                c->data[0], c->data[10], c->data[20], c->data[30],
                c->data[ndata-2], c->data[ndata-1]);
#endif
    return ndata;
}

FILE *get_classfile_descriptor(const char *file, CLASS_INFO *info)
{
    int len, rec_size, typeok = 0, mex, m, a;
    char *rec = NULL;
    FILE *fp;

    rec_size = 4 * info->reclen;
    rec = (char *)malloc(rec_size);
    if (!rec) {
        warning("Cannot allocate record of size %d.\n", rec_size);
        return NULL;
    }
#ifdef ALLOC
    Rprintf("%p -- allocated %d characters for rec\n", rec, rec_size);
#endif

    fp = fopen(file, "r");                  /* open the file for reading  */
    if (!fp) {                              /* test for failure           */
        warning("Cannot open file '%s'.\n", file);
        free(rec);
#ifdef ALLOC
        Rprintf("%p -- free\n", rec);
#endif
        return NULL;
    }

    len = fread(rec, sizeof(char), CLASS_BLK_SIZE, fp);
    if (len != CLASS_BLK_SIZE) {
        warning("Cannot read %d byte record from file '%s'.\n", CLASS_BLK_SIZE, file);
        fclose(fp);
        free(rec);
#ifdef ALLOC
        Rprintf("%p -- free\n", rec);
#endif
        return NULL;
    }

    if (info) {
        memcpy(info->fic,       &rec[ 0], 4); info->fic[3] = '\0';
        if (strncmp(info->fic, "1A", 2) == 0) {
            typeok = 1;
            memcpy(&(info->nextbl), &rec[ 4], 4);
            memcpy(&(info->nie),    &rec[ 8], 4);
            memcpy(&(info->nex),    &rec[12], 4);
            memcpy(&(info->first),  &rec[16], 4);
#ifdef DEBUG
            Rprintf("fic '%s'  (Type of data '1A' -> IEEE)\n", info->fic);
            Rprintf("nextbl '%d' (Next free block)\n", info->nextbl);
            Rprintf("nie '%d' (No of index entries)\n", info->nie);
            Rprintf("nex '%d' (No if index extensions)\n", info->nex);
            Rprintf("first '%d' (First block of 1st index extension)\n", info->first);
#endif
        }
        if (strncmp(info->fic, "2A", 2) == 0) {
            typeok = 2;
            memcpy(&(info->reclen), &rec[4], 4);
            memcpy(&(info->kind), &rec[8], 4);
            memcpy(&(info->vind), &rec[12], 4);
            memcpy(&(info->lind), &rec[16], 4);
            memcpy(&(info->flags), &rec[20], 4);

            memcpy(&(info->xnext), &rec[24], 8);
            memcpy(&(info->nextrec), &rec[32], 8);
            memcpy(&(info->nextword), &rec[40], 4);

            memcpy(&(info->lex1), &rec[44], 4);
            memcpy(&(info->nex), &rec[48], 4);
            memcpy(&(info->gex), &rec[52], 4);

            mex = info->nex;
            if (mex > 10) mex = 10;
            a = 56;
            for (m = 0; m < mex; m++) {
                memcpy(&(info->aex[m]), &rec[a], 8); a += 8;
            }

#ifdef DEBUG
            Rprintf("reclen '%d' (Record length [words])\n", info->reclen);
            Rprintf("kind '%d' (File kind)\n", info->kind);
            Rprintf("vind '%d' (Index version)\n", info->vind);
            Rprintf("lind '%d' (Index length [words])\n", info->lind);
            Rprintf("flags '%d' (Bit flags)\n", info->flags);
            Rprintf("xnext '%lld' (Next available entry number)\n", info->xnext);
            Rprintf("nextrec '%lld' (Next record which contains free space [recorrd])\n", info->nextrec);
            Rprintf("nextword '%d' (Next free word in this record [word])\n", info->nextword);
            Rprintf("lex1 '%d' (Length of first extension index [entries])\n", info->lex1);
            Rprintf("nex '%d' (No of index extensions)\n", info->nex);
            Rprintf("gex '%d' (Extension growth rule)\n", info->gex);
            for (m = 0; m < mex; m++) {
                Rprintf("  %d   aex=%lli record\n", m, info->aex[m]);
            }
#endif
            info->first = (int)info->lex1;
        }
        info->type = typeok;
    }
    if (info->nex != 1 && typeok == 1) {
        warning("Cannot read %d index extensions from file '%s' of type 1.", info->nex, file);
        fclose(fp);
        free(rec);
#ifdef ALLOC
        Rprintf("%p -- free\n", rec);
#endif
        return NULL;
    }
    free(rec);
#ifdef ALLOC
    Rprintf("%p -- free\n", rec);
#endif
    return fp;
}

int get_class_v1_data(FILE *fp, int curr_nbl, CLASS_INFO *info, SEXP a)
{
    int k, n, len, nbl, ns, nc;
    static char block[CLASS_BLK_SIZE], *sbl = NULL;
    const char *datetime = NULL;
    double restf, fres, LO, lam, bet, rchan; // vres
    bool spectrum;
    SEXP scan, classattrib;
    SEXP head, freq, data;
    SEXP nam;
    SEXP id, scanno, target, line, RA, Dec, f0, fLO, df, vs, dt, tsys, utc;

    CLASS_SECTION cobs;

    n = curr_nbl + 1;
    ns = 0;
    while (n < info->nextbl) {
        len = fread(block, sizeof(char), CLASS_BLK_SIZE, fp);
        nbl = check_block(st, info->first, n+1);
        if (nbl) {
            fill_class_v1_obs(&cobs, block);
#ifdef DEBUG
            Rprintf("%7d %12s %12s %d %d %d %d %d %d %d %d %d (%d %d %d %d | %d %d %d %d | %d %d %d %d)\n",
                    n,
                    st[nbl].xsourc, st[nbl].xline, st[nbl].xkind,
                    cobs.nbl, cobs.bytes, cobs.adr, cobs.nhead, cobs.len, cobs.ientry, cobs.nsec, cobs.obsnum,
                    cobs.sec_cod[0], cobs.sec_cod[1], cobs.sec_cod[2], cobs.sec_cod[3],
                    cobs.sec_adr[0], cobs.sec_adr[1], cobs.sec_adr[2], cobs.sec_adr[3],
                    cobs.sec_len[0], cobs.sec_len[1], cobs.sec_len[2], cobs.sec_len[3]);
#endif
            PROTECT(scan = allocVector(VECSXP, 3)); // head, freq, data
            classattrib = PROTECT(allocVector(STRSXP, 1));
            SET_STRING_ELT(classattrib, 0, mkChar("spectrum"));
            setAttrib(scan, R_ClassSymbol, classattrib);
            UNPROTECT(1);
            PROTECT(head = allocVector(VECSXP, 13));
            PROTECT(nam = allocVector(STRSXP, 13)); // names attribute (column names)
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
            UNPROTECT(1);

            PROTECT(id = allocVector(INTSXP, 1));         // 0
            PROTECT(scanno = allocVector(INTSXP, 1));     // 1
            PROTECT(target = allocVector(STRSXP, 1));     // 2
            PROTECT(line = allocVector(STRSXP, 1));       // 3
            PROTECT(RA = allocVector(REALSXP, 1));        // 4
            PROTECT(Dec = allocVector(REALSXP, 1));       // 5
            PROTECT(fLO = allocVector(REALSXP, 1));       // 6
            PROTECT(f0 = allocVector(REALSXP, 1));        // 7
            PROTECT(df = allocVector(REALSXP, 1));        // 8
            PROTECT(vs = allocVector(REALSXP, 1));        // 9
            PROTECT(dt = allocVector(REALSXP, 1));        // 10
            PROTECT(tsys = allocVector(REALSXP, 1));      // 11
            PROTECT(utc = allocVector(STRSXP, 1));        // 12

            INTEGER(id)[0] = n;
            INTEGER(scanno)[0] = st[nbl].xscan;
            trim(st[nbl].xsourc);
            SET_STRING_ELT(target, 0, mkChar(st[nbl].xsourc));
            trim(st[nbl].xline);
            SET_STRING_ELT(line, 0, mkChar(st[nbl].xline));
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
            UNPROTECT(13);

            sbl = (char *)malloc(cobs.nbl * sizeof(block));
            if (!sbl) {
                warning("CLASS: Not enough memory to allocate %d blocks.", cobs.nbl);
                fclose(fp);
                return -1;
            }
#ifdef ALLOC
            Rprintf("%p -- allocated %d blocks of size %ld\n", sbl, cobs.nbl, sizeof(block));
#endif
            /* copy the already read block into the sbl */
            memcpy(sbl, block, 512);
            if (cobs.nbl > 1) { /* and the remaining, if any */
                len = fread(&sbl[512], sizeof(char), (cobs.nbl - 1)*512, fp);
                if (len != (cobs.nbl - 1)*512) {
                    warning("CLASS: Read only %d (of %d) bytes.", len, (cobs.nbl - 1)*512);
                    free(sbl);
#ifdef ALLOC
                    Rprintf("%p -- free\n", sbl);
#endif
                    fclose(fp);
                    return -1;
                }
                n += cobs.nbl - 1;
            }
            /* Scan all nsec headers */
            for (k = 0; k < cobs.nsec; k++) {
                fill_class_header(cobs.sec_cod[k], cobs.sec_adr[k], cobs.sec_len[k],
                                  sbl, cobs.nbl*512, &st[cobs.obsnum]);
            }

            spectrum = (st[cobs.obsnum].xkind == 0);
            if (spectrum) {
              rchan = (st[cobs.obsnum]).s.rchan;
              restf = (st[cobs.obsnum]).s.restf;
              LO = ((st[cobs.obsnum]).s.restf + (st[cobs.obsnum]).s.image)/2.0;
              fres = (st[cobs.obsnum]).s.fres;
            } else {
              rchan = (st[cobs.obsnum]).u.rpoin;
              restf = (st[cobs.obsnum]).u.tref;
              LO = ((st[cobs.obsnum]).u.freq + (st[cobs.obsnum]).u.cimag)/2.0;
              fres = (st[cobs.obsnum]).u.tres;
              /*
              Rprintf("scan is continuum: %f %f %f %f\n", rchan, restf, LO, fres);
              Rprintf("scan is continuum: %f %f\n", (st[cobs.obsnum]).u.aref, 
                                                    (st[cobs.obsnum]).u.apos);
              */
            }
            lam = (st[cobs.obsnum]).p.lam;
            bet = (st[cobs.obsnum]).p.bet;
            lam += (st[cobs.obsnum]).p.lamof/cos(bet);
            bet += (st[cobs.obsnum]).p.betof;
            REAL(RA)[0] = lam*180.0/M_PI;
            REAL(Dec)[0] = bet*180.0/M_PI;
            REAL(fLO)[0] = LO;
            REAL(f0)[0] = restf;
            REAL(df)[0] = fres;
            REAL(vs)[0] = (st[cobs.obsnum]).s.voff;
            REAL(dt)[0] = (st[cobs.obsnum]).g.time;
            REAL(tsys)[0] = (st[cobs.obsnum]).g.tsys;
            datetime = obstime((st[cobs.obsnum]).xdobs + 60549, (st[cobs.obsnum]).g.ut);
            SET_STRING_ELT(utc, 0, mkChar(datetime));

            /* Scan data */
            nc = fill_class_v1_data(&cobs, sbl, cobs.nbl*512, &st[cobs.obsnum]);
            if (nc > 0) {
                PROTECT(freq = allocVector(REALSXP, nc));
                PROTECT(data = allocVector(REALSXP, nc));
                for (k = 0; k < nc; k++) {
                    if (fres != 0.0) {
                        REAL(freq)[k] = ((double)(k+1)-rchan)*fres + restf;
                    } else {
                        REAL(freq)[k] = (double)(k+1);
                    }
                    REAL(data)[k] = (st[cobs.obsnum]).data[k];
                }
                free((st[cobs.obsnum]).data);
#ifdef ALLOC
                Rprintf("%p -- free\n", (st[cobs.obsnum]).data);
#endif
            } else {
                PROTECT(freq = R_NilValue);
                PROTECT(data = R_NilValue);
            }
            PROTECT(nam = allocVector(STRSXP, 3));
            SET_STRING_ELT(nam, 0, mkChar("head"));
            SET_STRING_ELT(nam, 1, mkChar("freq"));
            SET_STRING_ELT(nam, 2, mkChar("data"));
            namesgets(scan, nam);
            UNPROTECT(1);

            SET_VECTOR_ELT(scan, 0, head);
            SET_VECTOR_ELT(scan, 1, freq);
            SET_VECTOR_ELT(scan, 2, data);
            UNPROTECT(3);

            SET_VECTOR_ELT(a, ns, scan);
            UNPROTECT(1);
            ns++;
#ifdef DEBUG
            Rprintf("Fill class data done for block %d.\n", n);
#endif
            free(sbl);
#ifdef ALLOC
            Rprintf("%p -- free\n", sbl);
#endif
        }
        n++;
    }

#ifdef DEBUG
    Rprintf("Return get_class_v1_data() at n = %d\n", n);
#endif
    return ns;
}

int get_class_v2_data(FILE *fp, int nscans, CLASS_INFO *info, int extno, SEXP a, int ns0)
{
    int k, n, pos, datapos, secpos, len, m, nmax, nbl, nc, ns;
    int section_size, code, length, datasize;
    char block[CLASS_BLK_SIZE];
    char *section, *databl;
    const char *datetime = NULL;
    double restf, fres, LO, lam, bet, rchan;
    SEXP scan, classattrib;
    SEXP head, freq, data;
    SEXP nam;
    SEXP id, scanno, target, line, RA, Dec, f0, fLO, df, vs, dt, tsys, utc;

    CLASS_SECTION_v2 cobs;
    CLASS *o;

    n = 0;
    ns = ns0;
    while (n < nscans) {
        o = &st[ns];
        pos = (o->xbloc - 1)*info->reclen + o->xword - 1;
#ifdef DEBUG
        Rprintf("n=%4d  Rec:%d Wor:%d  -> pos=%d words\n", n, o->xbloc, o->xword, pos);
#endif
        fseek(fp, 4*pos, SEEK_SET);
        len = fread(block, sizeof(char), CLASS_BLK_SIZE, fp);
        //        if (nbl) {
        fill_class_v2_obs(&cobs, block);
#ifdef DEBUG
        Rprintf("%7d %12s %12s %d %d (%d %d %d %d | %d %d %d %d | %d %d %d %d)\n",
                n,
                st[ns].xsourc, st[ns].xline, st[ns].xkind,
                cobs.nsec,
                cobs.sec_cod[0], cobs.sec_cod[1], cobs.sec_cod[2], cobs.sec_cod[3],
                cobs.sec_adr[0], cobs.sec_adr[1], cobs.sec_adr[2], cobs.sec_adr[3],
                cobs.sec_len[0], cobs.sec_len[1], cobs.sec_len[2], cobs.sec_len[3]);
#endif
        PROTECT(scan = allocVector(VECSXP, 3)); // head, freq, data
        classattrib = PROTECT(allocVector(STRSXP, 1));
        SET_STRING_ELT(classattrib, 0, mkChar("spectrum"));
        setAttrib(scan, R_ClassSymbol, classattrib);
        UNPROTECT(1);
        PROTECT(head = allocVector(VECSXP, 13));
        PROTECT(nam = allocVector(STRSXP, 13)); // names attribute (column names)
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
        UNPROTECT(1);

        PROTECT(id = allocVector(INTSXP, 1));         // 0
        PROTECT(scanno = allocVector(INTSXP, 1));     // 1
        PROTECT(target = allocVector(STRSXP, 1));     // 2
        PROTECT(line = allocVector(STRSXP, 1));       // 3
        PROTECT(RA = allocVector(REALSXP, 1));        // 4
        PROTECT(Dec = allocVector(REALSXP, 1));       // 5
        PROTECT(fLO = allocVector(REALSXP, 1));       // 6
        PROTECT(f0 = allocVector(REALSXP, 1));        // 7
        PROTECT(df = allocVector(REALSXP, 1));        // 8
        PROTECT(vs = allocVector(REALSXP, 1));        // 9
        PROTECT(dt = allocVector(REALSXP, 1));        // 10
        PROTECT(tsys = allocVector(REALSXP, 1));      // 11
        PROTECT(utc = allocVector(STRSXP, 1));        // 12

        INTEGER(id)[0] = n;
        INTEGER(scanno)[0] = st[ns].xscan;
        trim(st[ns].xsourc);
        SET_STRING_ELT(target, 0, mkChar(st[ns].xsourc));
        trim(st[ns].xline);
        SET_STRING_ELT(line, 0, mkChar(st[ns].xline));
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
        UNPROTECT(13);

        nmax = cobs.nsec;
        if (nmax > MAX_CLASS_SECT) nmax = MAX_CLASS_SECT;

        for (m = 0; m < nmax; m++) {
            section_size = 4*cobs.sec_len[m];
            section = (char *)malloc(section_size * sizeof(char));
            if (!section) {
                warning("Cannot allocate section (size=%d) for CLASS file.", section_size);
                fclose(fp);
                return -2;
            }
#ifdef ALLOC
            Rprintf("%p -- allocated %d characters for section\n", section, section_size);
#endif
            /* For now we read each section separately - to speed up all
               sections could be read in one go*/
            secpos = 4*(pos + cobs.sec_adr[m] - 1);
            fseek(fp, secpos*sizeof(char), SEEK_SET);
            len = fread(section, sizeof(char), section_size, fp);
            if (len != section_size) {
                warning("Cannot read %d word index from %d CLASS file.", len, section_size);
                free(section);
#ifdef ALLOC
                Rprintf("%p -- free\n", section);
#endif
                fclose(fp);
                return -2;
            }
            code = cobs.sec_cod[m];
            length = cobs.sec_len[m];
#ifdef DEBUG
            printf("Sect[%d]=%3d at adr %d of length %d.\n", m, code, secpos, section_size);
#endif
            fill_class_header(code, 1, length, section, section_size, o);
            free(section);
#ifdef ALLOC
            Rprintf("%p -- free\n", section);
#endif
        }

        rchan = (st[ns]).s.rchan;
        restf = (st[ns]).s.restf;
        LO = ((st[ns]).s.restf + (st[ns]).s.image)/2.0;
        fres = (st[ns]).s.fres;
        lam = (st[ns]).p.lam;
        bet = (st[ns]).p.bet;
        lam += (st[ns]).p.lamof/cos(bet);
        bet += (st[ns]).p.betof;
        REAL(RA)[0] = lam*180.0/M_PI;
        REAL(Dec)[0] = bet*180.0/M_PI;
        REAL(fLO)[0] = LO;
        REAL(f0)[0] = restf;
        REAL(df)[0] = fres;
        REAL(vs)[0] = (st[ns]).s.voff;
        REAL(dt)[0] = (st[ns]).g.time;
        REAL(tsys)[0] = (st[ns]).g.tsys;
        datetime = obstime((st[ns]).xdobs + 60549, (st[ns]).g.ut);
        SET_STRING_ELT(utc, 0, mkChar(datetime));

        /* Scan data */
        datasize = 4*cobs.ldata;
        datapos = 4*(pos + cobs.adata - 1);
        databl = (char *)malloc(datasize * sizeof(char));
        if (!databl) {
            warning("Cannot allocate data block (size=%d) for CLASS file.", datasize);
            fclose(fp);
            return -3;
        }
#ifdef ALLOC
        Rprintf("%p -- allocated %d characters for datatbl\n", databl, datasize);
#endif
        fseek(fp, datapos*sizeof(char), SEEK_SET);
        len = fread(databl, sizeof(char), datasize, fp);
        if (len != datasize) {
            warning("Read only %d of %d bytes data block CLASS file.", len, datasize);
            free(databl);
#ifdef ALLOC
            Rprintf("%p -- free\n", databl);
#endif
            fclose(fp);
            return -2;
        }
        nc = fill_class_v2_data(&cobs, databl, datasize, o);
        if (nc > 0) {
            PROTECT(freq = allocVector(REALSXP, nc));
            PROTECT(data = allocVector(REALSXP, nc));
            for (k = 0; k < nc; k++) {
                if (fres != 0.0) {
                    REAL(freq)[k] = ((double)(k+1)-rchan)*fres + restf;
                } else {
                    REAL(freq)[k] = (double)(k+1);
                }
                REAL(data)[k] = (st[ns]).data[k];
            }
            free((st[ns]).data);
#ifdef ALLOC
            Rprintf("%p -- free\n", (st[ns]).data);
#endif
        } else {
            PROTECT(freq = R_NilValue);
            PROTECT(data = R_NilValue);
        }
        PROTECT(nam = allocVector(STRSXP, 3));
        SET_STRING_ELT(nam, 0, mkChar("head"));
        SET_STRING_ELT(nam, 1, mkChar("freq"));
        SET_STRING_ELT(nam, 2, mkChar("data"));
        namesgets(scan, nam);
        UNPROTECT(1);

        SET_VECTOR_ELT(scan, 0, head);
        SET_VECTOR_ELT(scan, 1, freq);
        SET_VECTOR_ELT(scan, 2, data);
        UNPROTECT(3);

        SET_VECTOR_ELT(a, ns, scan);
        UNPROTECT(1);
        ns++;
#ifdef DEBUG
        Rprintf("Fill class data done for block %d.\n", n);
#endif
        free(databl);
#ifdef ALLOC
        Rprintf("%p -- free\n", databl);
#endif
        //        }
        n++;
    }
#ifdef DEBUG
    printf("Return get_class_v2_data() at n = %d\n", n);
#endif
    return ns;
}

void set_classfile_type(const char *file, CLASS_INFO *info)
{
    int len, record_length = 0;
    FILE *fp;
    static char msg[80];
    static char type[4], reclen[4];

    if (!info) return;

    info->type = 0;

    fp = fopen(file, "r");                  /* open the file for reading  */
    if (!fp) {                              /* test for failure           */
	sprintf(msg, "Cannot open file '%s'.", file);
	error(msg);
    }

    len = fread(type, sizeof(char), 4, fp);
    if (len != 4) {
        warning("Cannot read %d byte type from file '%s'.\n", 4, file);
        fclose(fp);
        return;
    }

    if (strncmp(type, "1A", 2) == 0) {
        info->type = 1;
        info->reclen = 128;
    } else if (strncmp(type, "2A", 2) == 0) {
        info->type = 2;
        len = fread(reclen, sizeof(char), 4, fp);
        if (len != 4) {
            warning("Cannot read %d byte reclen from file '%s'.\n", 4, file);
            fclose(fp);
            return;
        }
        memcpy(&record_length, reclen, 4);
        info->reclen = record_length;
    } else {
        fclose(fp);
        return;
    }

    fclose(fp);

#ifdef DEBUG
    printf("Found CLASS file of type='%s' and record length=%d words.\n", type, record_length);
#endif

    return;
}

//' Read a GILDAS/CLASS single dish data file
//'
//' Given a filename, open the file and scan it for single dish spectra
//' or continuum scans. Each encountered scan is returned as a list with
//' a head, freq and data section. For continuum scans, the frequency
//' vector will simply be a running index. All the individual lists are
//' combined into one major list, which is returned.
//' @param filename name of the GILDAS file including path to be opened
//' @return list of length n, where n is the number of scans found
// [[Rcpp::export]]
SEXP readClass(SEXP filename)
{
    const char *s;
    int ierr, n, nbl, nspec, nex, m, ns0;
    FILE *fp;
    CLASS_INFO info;
    SEXP classattrib, ret;

    s = CHAR(STRING_ELT(filename, 0));
    Rprintf("CLASS filename = %s\n", s);

    info.type = 0;
    set_classfile_type(s, &info);
    if (info.type == 0) {
        warning("Unknown CLASS file type in %s.", s);
        return 0;
    }

#ifdef DEBUG
    printf("File: '%s'  type=%d\n", s, info.type);
#endif

    fp = get_classfile_descriptor(s, &info);
    if (fp) {
        if (info.type == 1) {
            nspec = get_class_v1_file_listing(fp, &info);
            if (nspec <= 0) {
                if (st) free(st);
#ifdef ALLOC
                Rprintf("%p -- free\n", st);
#endif
                st = NULL;
                nst = 0;
                fclose(fp);
                return R_NilValue;
            }
            n = nspec/4+2; // number of 512 byte block
            PROTECT(ret = allocVector(VECSXP, nspec)); // a list with nspec elements
            ierr = get_class_v1_data(fp, n, &info, ret);
            classattrib = PROTECT(allocVector(STRSXP, 1));
            SET_STRING_ELT(classattrib, 0, mkChar("spectra"));
            setAttrib(ret, R_ClassSymbol, classattrib);
            UNPROTECT(2);
            fclose(fp);
            if (st != NULL) free(st);
#ifdef ALLOC
            Rprintf("%p -- free\n", st);
#endif
            return ret;
        } else if (info.type == 2) {
            nex = info.nex;
            nspec = 0;
            for (m = 0; m < nex; m++) {
                nspec += get_class_v2_file_listing(fp, &info, m, 0);
            }
            PROTECT(ret = allocVector(VECSXP, nspec)); // a list with nspec elements
            ns0 = 0;
            for (m = 0; m < nex; m++) {
                nspec = get_class_v2_file_listing(fp, &info, m, 1);
                if (nspec <= 0) {
                    if (st) free(st);
#ifdef ALLOC
                    Rprintf("%p -- free\n", st);
#endif
                    st = NULL;
                    nst = 0;
                    fclose(fp);
                    return R_NilValue;
                }
                ierr = get_class_v2_data(fp, nspec, &info, m, ret, ns0);
                ns0 += nspec;
                if (st != NULL) free(st);
#ifdef ALLOC
                Rprintf("%p -- free\n", st);
#endif
            }
            classattrib = PROTECT(allocVector(STRSXP, 1));
            SET_STRING_ELT(classattrib, 0, mkChar("spectra"));
            setAttrib(ret, R_ClassSymbol, classattrib);
            UNPROTECT(2);
            fclose(fp);
            return ret;
        }
    } else {
        return R_NilValue;
    }
}
