options <- new.env()

.onLoad <- function(libname, pkgname) {
    # packageStartupMessage('This is Rdrp\n', domain = NULL, appendLF = TRUE)
    # print(ls(envir=options))
    options$system <- "frequency"
    options$position.tolerance <- 1/3600  # 1 arcsec
    options$frequency.tolerance <- 1.0e6  # 1 MHz
}
