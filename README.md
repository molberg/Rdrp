# Rdrp

The **Rdrp** package provides commands for the reduction of spectral line
data originating from single dish radio astronomical observations. In
that sense it is similar to programs like CLASS, which is part of the
[GILDAS](https://www.iram.fr/IRAMFR/GILDAS/) system. But
contrary to CLASS (and other alternatives), **Rdrp** is built as a
package for the **R programming language**, thus building on a huge
reservoir of other numerical and statistical tools. Also, the user
gets the high quality, powerful built-in graphics from R for free.

The principle design is as follows: data are organized in memory in
the form of lists (of class *spectrum*), which have three componnents:

 * *head* header information for the spectrum, itself a list
 * *freq* a numeric vector giving the frequencies of all spectral
   channels
 * *data* a numeric vector giving the data content of all spectral channels

The available routines to read data from a FITS file or a CLASS
formatted binary file, typically return lists of objects, each of
which will be of the format described above. Helper functions
*getHead*, *getFreq* and *getData* are available to return all headers as a data
frame (with one row per spectrum) and all the frequency (or data)
vectors as a numeric matrix, where spectral channels run along rows
and the columns correspond to an individual spectrum. Such a list of
lists would be of class *spectra*, i.e. plural of *spectrum*.

Here is an example of what a short session may look like:

``` r
    library(Rdrp)
	# 
    # Here is how you would construct a very simple, fake spectrum:
    head <- list(target="my target", ra=1.0, dec=2.0, f0=1421.0)
    freq <- head$f0 + seq(-100,100)*0.1  # +-10 MHz around centre frequency
    data <- 0.1*rnorm(length(freq)) + exp(-((freq-head$f0)/10)^2)
    S <- list(head=head, freq=freq, data=data)
    class(S) <- "spectrum"
    plot(S)             # this will use plot.spectrum(...)
    #
    # a somewhat more realistic session, using data from APEX
    assign("system","velocity", Rdrp::options) # work in velocity space
    L <- readClass("mydata.apex")              # class(L) is 'spectra'
    H <- getHead(L)
    print(H)            # take a look at the header information
    i <- which(H$target == "IC348" & H$line == "CO(3-2)")
    L <- L[i]           # only keep spectra with given target, line
    S <- L[[1]]         # get the first spectrum, class(S) is 'spectrum'
    plot(S)             # and plot it
    A <- average(L)     # form the average
    # we expect a spectral line between -20 and +20 km/s
    linemask <- mask(A, c(-20,20))
    # fit a second order baseline
    bl <- baseline(A, order=2, mask=linemask)
    plot(A, type='s')                 # plot the average spectrum ...
    lines(velocity(A), bl, col='red') # ... and fitted baseline
    A$data <- A$data - bl             # subtract baseline
    plot(A, type='s', xlab="velocity [km/s]", ylab=expression(T[A]))
```
