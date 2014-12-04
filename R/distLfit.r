# fit many different distributions via linear moments and return goodness of fit values
# Berry Boessenkool, Sept 2014
# Main function calling   distLgof   and   distLgofplot or distLplot

distLfit <- function(
dat,          # Vector with values
datname,      # Character string for main, xlab etc
speed=TRUE,   # If TRUE, several distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}
ks=TRUE,      # Include ks.test results in dlf$gof? Computing is much faster when FALSE
selection,    # Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}. Overrides speed.
gofProp=1,    # Upper proportion of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
gofComp=FALSE,# If TRUE, plots a comparison of the ranks of different GOF-methods and sets plot to FALSE
progbars=TRUE,# Show progress bars for each loop?
time=TRUE,    # \code{\link{message}} execution time?
plot=TRUE,    # Should a histogram with densities be plotted?
cdf=FALSE,    # If TRUE, plot cumulated DF instead of probability density
legargs=NULL, # List of arguments passed to \code{\link{legend}} except for legend and col
histargs=NULL,# List of arguments passed to \code{\link{hist}} except for x, breaks, col, xlim, freq
quiet=FALSE,  # Should \code{\link{rmse}} warn about NA removal?
... )         # Further arguments passed to \code{\link{distLplot}} if they are accepted there, else passed to \code{\link{lines}}, like lty, type, pch, ...
{
# preparation ------------------------------------------------------------------
StartTime <- Sys.time()
if(missing(datname)) datname <- deparse(substitute(dat))
if(length(gofProp)>1 | any(gofProp<0) | any(gofProp>1) ) stop("gofProp must be a single value between 0 and 1")
# Progress bars
if( require(pbapply,quietly=TRUE) & progbars ) lapply <- pbapply::pblapply
# checks:
if( ! is.numeric(dat) ) stop("dat must be numeric.")
if( length(gofProp) != 1  |  gofProp > 1  |  gofProp < 0 )
   stop("gofProp must be a single value between 0 and 1.")
# remove NAs, convert to vector:
dat <- as.numeric( dat[!is.na(dat)]  )
# possible distributions:
dn <- dist.list()
# Selection:
if( ! missing(selection) )
  {
  if(is.numeric(selection)) if(any(abs(selection)>length(dn)))
     stop("'selection' cannot be larger than", length(dn))
  names(dn) <- dn # so that selection can also be character string
  dn <- dn[selection]
  names(dn) <- NULL
  }
else
# remove some to save time and errors, see ?dist.list # gld added
if(speed) dn <- dn[ ! dn %in%
   c("aep4","cau","emu","gep","gld","kmu","kur","lmrq","sla","st3","texp")]
#
# Fit distribution parameters --------------------------------------------------
# L-Moments of sample  # package lmomco
mom <- lmoms(dat, nmom=5)
# estimate parameters for each distribution:    # this takes time!
if(progbars) print("Parameter estimation from linear moments:")
parameter <- lapply(dn, function(d) lmom2par(mom, type=d) )
# error catching:
if( length(parameter) != length(dn))
  {
  on.exit(warning("Some distributions could not be fitted. Possibly cau."))
  names(parameter) <- sapply(parameter, "[[", "type")
  }
else names(parameter) <- dn
# Goodness of Fit, output list, plot -------------------------------------------
output <- distLgof(list(dat=dat, datname=datname, gofProp=gofProp, parameter=parameter),
                   plot=FALSE, progbars=progbars, ks=ks, quiet=quiet)
# compare GOF
if(gofComp)
  {
  distLgofplot(output)
  plot <- FALSE
  }
if(plot) distLplot(dlf=output, cdf=cdf, legargs=legargs, histargs=histargs, ... )
if(time) message("execution took ", signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.\n")
return(invisible(output))
} # end of function
