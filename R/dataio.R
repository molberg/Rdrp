## A routine to look-up a variable in a FITS header
lookup <- function(hdr, keyword) {
    i <- which(hdr$key == keyword)
    if (length(i) == 0) {
        value <- NA
    } else {
        value <- hdr$value[i[1]]
    }
    value
}

#' Read FITS files from the Onsala 20m telescope.
#'
#' Take a FITS file for a single spectrum from the OSO 20m and return a spectrum
#' (i.e. a list consisting of header, frequency and data vectors).
#' @param fitsfiles file name(s) including path of FITS file(s) to convert
#' @return a list of lists with components head, freq and data
#' @export
readOSO20m <- function(fitsfiles) {
    L <- lapply(fitsfiles,
                function(filename) {
                    MHz <- 1.0e6
                    f <- FITSio::readFITS(file=filename, hdu = 1, phdu = 1)
                    even <- seq(2, length(f$hdr), by=2)
                    odd <- even-1
                    hdr <- data.frame(key=f$hdr[odd], value=f$hdr[even], stringsAsFactors=FALSE)
                    data <- as.numeric(f$imDat)
                    f0 <- f$axDat$crval[1]/MHz
                    f.rest <- as.numeric(lookup(hdr, "RESTFREQ"))/MHz
                    if (f0 == 0.0 & is.numeric(f.rest)) {
                        f0 <- f.rest
                    }
                    RA <- f$axDat$crval[2]
                    Dec <- f$axDat$crval[3]
                    offx <- f$axDat$cdelt[2]
                    offy <- f$axDat$cdelt[3]
                    freq <- f0+f$axDat$cdelt[1]*(seq(f$axDat$len[1])-f$axDat$crpix[1])/MHz
                    id <- as.integer(lookup(hdr, "SCAN-NUM"))
                    dt <- as.double(lookup(hdr, "OBSTIME"))
                    target <- lookup(hdr, "OBJECT")
                    #f1 <- as.numeric(as.double(lookup(hdr, "IMAGFREQ")))
                    df <- as.double(lookup(hdr, "CDELT1"))/MHz
                    if (df < 0.0) {
                        df <- -df
                        freq <- rev(freq)
                        data <- rev(data)
                    }
                    vs <- as.double(lookup(hdr, "VELO-LSR"))/1.0e3
                    T.sys <- as.double(lookup(hdr, "TSYS"))
                    T.rec <- as.double(lookup(hdr, "TREC"))
                    T.amb <- as.double(lookup(hdr, "TAMB"))
                    p.amb <- as.double(lookup(hdr, "PRESSURE"))/100.0
                    rel.hum <- as.double(lookup(hdr, "HUMIDITY"))
                    if (is.na(vs)) vs <- as.numeric(lookup(hdr, "VLSR"))/1.0e3
                    line <- lookup(hdr, "LINE")
                    date <- lookup(hdr, "DATE-OBS")
                    eta.mb <- as.double(lookup(hdr, "BEAMEFF"))
                    if (nchar(date) == 10) {
                        date <- paste(date, lookup(hdr, "UTC"), sep="T")
                    }
                    # tstamp <- as.POSIXlt(date, format="%Y-%m-%dT%H:%M:%S")
                    Az <- as.double(lookup(hdr, "AZIMUTH"))
                    El <- as.double(lookup(hdr, "ELEVATIO"))
                    head <- list(id=id, target=target, line=line,
                                 RA=RA, Dec=Dec, offx=offx, offy=offy, Az=Az, El=El,
                                 f0=f0, v.LSR=vs, eta.mb=eta.mb,
                                 T.amb=T.amb, p.amb=p.amb, rel.hum=rel.hum,
                                 T.sys=T.sys, df=df, dt=dt, observed.date=date)
                    sd <- list(head=head, freq=freq, data=data)
                    class(sd) <- "spectrum"
                    sd
                })
    class(L) <- "spectra"
    L
}

#' Read a FITS file from the Onsala SALSA student telescope.
#'
#' Take a FITS file for a single spectrum from SALSA and return a spectrum
#' (i.e. a list consisting of header, frequency and data vectors).
#' @param fitsfiles file name(s) including path of FITS file(s) to convert
#' @return a list of lists with components head, freq and data
#' @export
readSALSA <- function(fitsfiles) {
    L <- lapply(fitsfiles,
                function(filename) {
                    MHz <- 1.0e6
                    id <- as.integer(gsub("[^0-9]", "", basename(filename)))
                    f <- FITSio::readFITS(file=filename, hdu = 1, phdu = 1)
                    even <- seq(2, length(f$hdr), by=2)
                    odd <- even-1
                    hdr <- data.frame(key=f$hdr[odd], value=f$hdr[even], stringsAsFactors=FALSE)
                    freq <- (f$axDat$crval[1]+f$axDat$cdelt[1]*(seq(f$axDat$len[1])-f$axDat$crpix[1]))/MHz
                    bzero <- as.double(lookup(hdr, "BZERO"))
                    if (is.na(bzero)) bzero <- 0.0
                    bscale <- as.double(lookup(hdr, "BSCALE"))
                    if (is.na(bscale)) bscale <- 1.0
                    data <- as.numeric(f$imDat*bscale + bzero)
                    target <- lookup(hdr, "OBJECT")
                    if (is.na(target)) target = "unknown"
                    onx <- as.double(lookup(hdr, "CRVAL2"))
                    ony <- as.double(lookup(hdr, "CRVAL3"))
                    dt <- as.double(lookup(hdr, "OBSTIME"))
                    f1 <- as.double(lookup(hdr, "CRVAL1"))/MHz
                    f0 <- as.double(lookup(hdr, "RESTFREQ"))/MHz
                    if (is.na(f0)) f0 <- 1420.40575177
                    df <- as.double(lookup(hdr, "CDELT1"))/MHz
                    vs <- as.double(lookup(hdr, "VLSR"))
                    if (is.na(vs)) vs <- as.double(lookup(hdr, "VELO-LSR"))
                    tsys <- as.double(lookup(hdr, "TSYS"))
                    date <- lookup(hdr, "DATE-OBS")
                    # tstamp <- as.POSIXlt(date, format="%Y-%m-%dT%H:%M:%S")
                    head <- list(id=id, target=target, line="21cm",
                                 LII=onx, BII=ony,
                                 f0=f0, f1=f1, v.LSR=vs,
                                 dt=dt, observed.date=date)
                    sd <- list(head=head, freq=freq, data=data)
                    class(sd) <- "spectrum"
                    sd
                })
    class(L) <- "spectra"
    L
}
