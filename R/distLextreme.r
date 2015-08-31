# Extreme value statistics for flood risk estimation
# Berry Boessenkool, 2012 (first draft) - 2014 (main update)

# Input: vector with annual discharge maxima
# Ouput: discharge estimates for given return periods, parameter of distributions, quality of fits, plot with linear axis
# See German RclickHandbuch.wordpress.com Kapitel 15

distLextreme <- function(
dat,          # Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.
datname,      # Character string for main and output list.
dlf,          # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat!
selection=NULL, # Selection of distributions. Character vector with types as in \code{\link[lmomco]{lmom2par}}. Overrides speed.
RPs=c(2,5,10,20,50), # ReturnPeriods for which discharge is estimated
progbars=length(dat)>200, # Show progress bars for each loop? DEFAULT: TRUE if n>200
time=TRUE,    # \code{\link{message}} execution time?
plot=TRUE,    # Should the return periods and nbest fitted distributions be plotted by a call to \code{\link{distLextremePlot}}?
quiet=FALSE,  # Suppress notes?
... )         # Further arguments passed to \code{\link{distLextremePlot}} like order, lty, lwd, ...
{
StartTime <- Sys.time()
if(quiet) progbars <- FALSE
# input checks -----------------------------------------------------------------
if(missing(datname)) datname <- deparse(substitute(dat))
if(!missing(dlf)) if(!missing(dat)) if(dlf$dat != dat & !quiet)
   message("note in distLextreme: 'dat' differs from 'dlf$dat'. 'dat' is ignored.")
#
# fit distributions and calculate goodness of fits -----------------------------
if( missing(dlf) )  dlf <- distLfit(dat=dat, datname=datname, plot=FALSE, 
                selection=selection, time=FALSE, progbars=progbars, quiet=quiet)
#
# plot -------------------------------------------------------------------------
if(plot) dlf <- distLextremePlot(dlf=dlf, selection=selection, quiet=quiet, ...)
# output (discharge) values at return periods ----------------------------------
dn <- rownames(dlf$gof) # distribution names
if(progbars) sapply <- pbapply::pbsapply
if(progbars) message("Calculating return levels for return periods:")
returnlev <- sapply(dn, function(d) qlmomco(1-1/RPs, dlf$parameter[[d]])) 
# if length(RPs)==1, returnlev is only a vector,
if(is.null(dim(returnlev)))
     returnlev <- as.matrix(returnlev)        # so convert it to a matrix,
else returnlev <- t(returnlev)                # or else transpose it
# column names:
colnames(returnlev) <- paste0("RP.", RPs)
# add weighted estimate (by goodness of fit):
# todo...
# Add to output:
dlf$returnlev <- as.data.frame(returnlev)
if(time & !quiet) message("distLextreme execution took ", signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(dlf))
} # end of function
