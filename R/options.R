.onLoad <- function(libname, pkgname) {
    #drp.options <- new.env()
    #drp.options$system <- "velocity"
    #drp.options$position.tolerance <- 1/3600  # 1 arcsec
    #drp.options$frequency.tolerance <- 1.0e6  # 1 MHz
    drp.options = list(system="velocity")

    packageStartupMessage('This is Rdrp\n', domain = NULL, appendLF = TRUE)
    print(drp.options)
}
