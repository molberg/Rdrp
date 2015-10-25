#' Return headers of spectra.
#'
#' Turn head component of list of spectra into data.frame
#' @param L list of spectra
#' @return a dataframe containing header info for each spectrum, one per row.
#' @export
header <- function(L) {
    head <- lapply(L,
                   function(S) {
                       data.frame(S$head, stringsAsFactors = FALSE)
                   })
    do.call(rbind, head)
}

frequency <- function(L) {
    freq <- lapply(L,
                   function(S) {
                       as.numeric(S$freq)
                   })
    mat <- do.call(cbind, freq)
    as.matrix(mat)
}

data <- function(L) {
    data <- lapply(L,
                   function(S) {
                       as.numeric(S$data)
                   })
    mat <- do.call(cbind, data)
    as.matrix(mat)
}

#' Make a dataset from a list of individual spectra.
#'
#' Given a list of spectra (i.e. lists with head, freq and data components),
#' return a dataset (a list) which has a data.frame describing header information,
#' and two matrices (columns of frequency and data vectors).
#' @param L a list of lists with components head, freq, data
#' @return a list with components head (data.frame), freq (matrix) and data (matrix)
#' @export
makeDataset <- function(L) {
    h <- header(L)
    h$target <- as.factor(h$target)
    h$line <- as.factor(h$line)
    f <- frequency(L)
    d <- data(L)
    sd = list(head=h, freq=f, data=d)
    class(sd) <- "spectra"
    sd
}

#' Pick a subset of spectra.
#'
#' Given an integer vector of indices, construct a subset of a dataset by selecting
#' only the rows (columns) from the header (data) of the original dataset.
#' @param ds the original dataset
#' @param index the rows (columns) to pick
#' @return an object of the same kind as ds, but only keeping given rows (columns)
#' @export
pick <- function(ds, index, drop.levels=FALSE) {
    h <- ds$head[index,]
    if (drop.levels) droplevels(h)
    f <- ds$freq[,index]
    d <- ds$data[,index]
    sd <- list(head=h, freq=f, data=d)
    class(sd) <- "spectra"
    sd
}
