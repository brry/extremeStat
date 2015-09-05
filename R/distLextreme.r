# Extreme value statistics for flood risk estimation
# Berry Boessenkool, 2012 (first draft) - 2014 & 2015 (main updates)
# berry-b@gmx.de

# Input: vector with annual discharge maxima
# Ouput: discharge estimates for given return periods, parameter of distributions, quality of fits, plot with linear axis
# See German RclickHandbuch.wordpress.com Kapitel 15

distLextreme <- function(
dat,          # Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.
dlf=NULL,     # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat!
selection=NULL, # Selection of distributions. Character vector with types as in \code{\link[lmomco]{lmom2par}}.
truncate=0,     # Proportion truncated for censored quantile, see \code{\link{distLquantile}}.
RPs=c(2,5,10,20,50), # ReturnPeriods for which discharge is estimated
progbars=length(dlf$dat)>200, # Show progress bars for each loop? DEFAULT: TRUE if n>200
time=TRUE,    # \code{\link{message}} execution time?
plot=TRUE,    # Should the return periods and nbest fitted distributions be plotted by a call to \code{\link{distLextremePlot}}?
quiet=FALSE,  # Suppress notes?
... )         # Further arguments passed to \code{\link{distLextremePlot}} like order, lty, lwd, ...
{
StartTime <- Sys.time()
if(quiet) progbars <- FALSE
# fit distributions and calculate goodness of fits -----------------------------
if( is.null(dlf) )  dlf <- distLfit(dat=dat, datname=deparse(substitute(dat)), 
      plot=FALSE, selection=selection, time=FALSE, progbars=progbars, quiet=quiet)
# Equality check
if(!missing(dat) & !is.null("dlf")) if(any(dlf$dat != dat) & !quiet)
  on.exit(message("Note in distLextreme: 'dat' differs from 'dlf$dat'. 'dat' is ignored."))
#
# plot -------------------------------------------------------------------------
if(plot) dlf <- distLextremePlot(dlf=dlf, selection=selection, quiet=quiet, ...)
#
# output (discharge) values at return periods ----------------------------------
returnlev <- distLquantile(dlf=dlf, selection=selection, truncate=truncate, 
                           probs=1-1/RPs, empirical=FALSE, weighted=TRUE, trans=TRUE)
# column names:
colnames(returnlev) <- paste0("RP.", RPs)
# Add to output:
dlf$returnlev <- as.data.frame(returnlev)
if(time & !quiet) on.exit(message("distLextreme execution took ", 
  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds."))
return(invisible(dlf))
} # end of function
