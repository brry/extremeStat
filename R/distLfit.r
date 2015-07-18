# fit many different distributions via linear moments and return goodness of fit values
# Berry Boessenkool, Sept 2014, July 2015
# Main function calling   distLgof   and   distLgofPlot or distLplot

distLfit <- function(
dat,          # Vector with values
datname,      # Character string for main, xlab etc. DEFAULT: internal \code{substitute(dat)}
speed=TRUE,   # If TRUE, several distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}
ks=FALSE,     # Include ks.test results in \code{dlf$gof}? Computing is much faster when FALSE
selection=NULL, # Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}. Overrides speed.
gofProp=1,    # Upper proportion (0:1) of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
truncate=0,   # Number between 0 and 1. Censored \code{\link{distLquantile}}: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly.
gofComp=FALSE,# If TRUE, plots a comparison of the ranks of different GOF-methods and sets plot to FALSE
progbars=length(dat)>200,# Show progress bars for each loop? DEFAULT: TRUE if n > 200
time=TRUE,    # \code{\link{message}} execution time?
plot=TRUE,    # Should a histogram with densities be plotted?
cdf=FALSE,    # If TRUE, plot cumulated DF instead of probability density
legargs=NULL, # List of arguments passed to \code{\link{legend}} except for legend and col
histargs=NULL,# List of arguments passed to \code{\link{hist}} except for x, breaks, col, xlim, freq
quiet=FALSE,  # Suppress notes?
... )         # Further arguments passed to \code{\link{distLplot}} if they are accepted there, else passed to \code{\link{lines}}, like lty, type, pch, ...
{
# preparation ------------------------------------------------------------------
StartTime <- Sys.time()
if(missing(datname)) datname <- deparse(substitute(dat))
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
# checks:
if( ! is.numeric(dat) ) stop("dat must be numeric.")
if(!is.vector(dat) & !quiet) message("note in distLfit: dat was not a vector.")
if(length(gofProp)>1 | any(gofProp<0) | any(gofProp>1) )
  stop("gofProp must be a single value between 0 and 1.")
# remove NAs, convert to vector:
dat <- dat_full <- as.numeric( dat[!is.na(dat)]  )
# truncate (fit values only to upper part of values):
if(truncate!=0) dat <- sort(dat)[ -1:-(truncate*length(dat)) ]
# Check remaining sample size
if(length(dat) < 5) {message("note in distLfit: sample size (",length(dat),") is too small to fit parameters (<5).")
                     error_out <- as.list(selection) # this is very useful for distLquantile
                     names(error_out) <- selection  # since it keeps the columns if a selection is given
                     return(list(dat=dat, parameter=error_out, gof=NA))}
#
# possible distributions: ------------------------------------------------------
dn <- dist.list()
# Selection:
if( ! is.null(selection) )
  {
  if(is.numeric(selection)) if(any(abs(selection)>length(dn)))
     stop("'selection' cannot be larger than", length(dn), ".")
  names(dn) <- dn # so that selection can also be character string
  seldn <- !selection %in% dn
  if(any(seldn))
   {
   if(!quiet) message("note in distLfit: selection '", paste(selection[seldn], collapse="', '"),
   "' not available in lmomco::dist.list(), thus removed.")
   selection <- selection[!seldn]
   }
  dn <- dn[selection]
  names(dn) <- NULL
  }
else
# remove some to save time and errors, see ?dist.list # gld, gov and tri added
if(speed) dn <- dn[ ! dn %in%
   c("aep4","cau","emu","gep","gld","gov","kmu","kur","lmrq","sla","st3","texp","tri")]
#
# Fit distribution parameters --------------------------------------------------
# L-Moments of sample  # package lmomco
mom <- lmoms(dat, nmom=5)
# estimate parameters for each distribution:    # this takes time!
if(progbars) message("Parameter estimation from linear moments:")
parameter <- lapply(dn, function(d) lmom2par(mom, type=d) )
# error catching:
if( length(parameter) != length(dn))
  {
  if(!quiet) on.exit(message("note in distLfit: Some distributions could not be fitted. Possibly cau."))
  names(parameter) <- sapply(parameter, "[[", "type")
  }
else names(parameter) <- dn
#
# Goodness of Fit, output list, plot -------------------------------------------
output <- distLgof(list(dat=dat, datname=datname, gofProp=gofProp, parameter=parameter),
                   plot=FALSE, progbars=progbars, ks=ks, quiet=quiet)
# compare GOF
if(gofComp)
  {
  distLgofPlot(output, quiet=quiet)
  plot <- FALSE
  }
if(plot) output <- distLplot(dlf=output, cdf=cdf, legargs=legargs, histargs=histargs, ... )
if(!plot) output$coldist <- rainbow2(if(is.null(selection)) 5 else length(selection))
# truncation value
output$truncate <- truncate
output$dat_full <- dat_full # non-truncated data
if(time & !quiet) message("distLfit execution took ", signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(output))
} # end of function
