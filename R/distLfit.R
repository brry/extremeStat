#' Fit distributions via linear moments
#' 
#' Fit several distributions via linear moments, plot histogram and
#' distribution densities \emph{or} ecdf with cumulated probability.
#' Also returns goodness of fit values.
#' This is the main fitting function calling   distLgof   and   distLgofPlot or distLplot
#'
#' @details Fits parameters via \code{\link[lmomco]{lmom2par}} in the package \code{lmomco}
#' @return List as explained in \code{\link{extremeStat}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014 + July 2015
#' @seealso \code{\link{distLgof}}, \code{\link{distLplot}}.
#'          \code{\link[extRemes]{fevd}} in the package \code{extRemes},
#'          \code{\link[MASS]{fitdistr}} in the package \code{MASS}.
#' @keywords hplot dplot distribution
#' @export
#' @importFrom lmomco dist.list lmoms lmom2par
#' @importFrom berryFunctions rainbow2
#' 
#' @examples
#' 
#' data(annMax)
#' # basic usage on real data (annual discharge maxima in Austria)
#' dlf <- distLfit(annMax)
#' str(dlf, max.lev=2)
#' distLprint(dlf)
#' 
#' # arguments that can be passed:
#' distLfit(annMax, lty=2, col=3, legargs=list(lwd=3), main="booh!")
#' set.seed(42)
#' dlf_b <- distLfit(rbeta(100, 5, 2), nbest=10, legargs=c(x="left"))
#' distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"))
#' distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"), order=TRUE)
#' distLplot(dlf_b, coldist=c("orange",3:6), lty=1:3) # lty is recycled
#' distLplot(dlf_b, cdf=TRUE)
#' distLplot(dlf_b, cdf=TRUE, histargs=list(do.points=FALSE), sel="nor")
#' 
#' 
#' # Goodness of Fit is computed by RMSE, see first example of ?distLgof
#' 
#' # logarithmic axes:
#' set.seed(1)
#' y <- 10^rnorm(100, mean=2, sd=0.3) # if you use 1e4, distLgof will be much slower
#' hist(y, breaks=20)
#' berryFunctions::logHist(y, col=8)
#' dlf <- distLfit(log10(y), breaks=50)
#' distLplot(dlf, breaks=50, log=TRUE)
#' 
#' \dontrun{
#' # this takes a while, as it tries to fit all 30 distributions:
#' d_all <- distLfit(annMax, gofProp=1, speed=FALSE, plot=FALSE) # 35 sec
#' distLprint(d_all)
#' distLplot(d_all, nbest=22, coldist=grey(1:22/29), xlim=c(20,140))
#' distLplot(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
#' d_all$gof
#' }
#' 
#' @param dat Vector with values
#' @param datname Character string for main, xlab etc. DEFAULT: internal \code{substitute(dat)}
#' @param speed If TRUE, several distributions are omitted, for the reasons shown in \code{lmomco::\link[lmomco]{dist.list}()}. DEFAULT: TRUE
#' @param ks Include ks.test results in \code{dlf$gof}? Computing is much faster when FALSE. DEFAULT: FALSE
#' @param selection Selection of distributions. Character vector with types as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL
#' @param gofProp Upper proportion (0:1) of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail. DEFAULT: 1
#' @param weightc Named custom weights for each distribution, see \code{\link{distLgof}}. DEFAULT: NA
#' @param truncate Number between 0 and 1. POT Censored \code{\link{distLquantile}}: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly. DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, set this only via truncate, see Details of \code{\link{q_gpd}}. DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param gofComp If TRUE, plots a comparison of the ranks of different GOF-methods and sets plot to FALSE. DEFAULT: FALSE
#' @param progbars Show progress bars for each loop? DEFAULT: TRUE if n > 200
#' @param time \code{\link{message}} execution time? DEFAULT: TRUE
#' @param plot Should a histogram with densities be plotted? DEFAULT: TRUE
#' @param cdf If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE
#' @param legargs List of arguments passed to \code{\link{legend}} except for legend and col. DEFAULT: NULL
#' @param histargs List of arguments passed to \code{\link{hist}} except for x, breaks, col, xlim, freq. DEFAULT: NULL
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param ssquiet Suppress sample size notes? DEFAULT: quiet
#' @param \dots Further arguments passed to \code{\link{distLplot}} if they are accepted there, else passed to \code{\link{lines}}, like lty, type, pch, ...
#' 
distLfit <- function(
dat,
datname,
speed=TRUE,
ks=FALSE,
selection=NULL,
gofProp=1,
weightc=NA,
truncate=0,
threshold=berryFunctions::quantileMean(dat, truncate),
gofComp=FALSE,
progbars=length(dat)>200,
time=TRUE,
plot=TRUE,
cdf=FALSE,
legargs=NULL,
histargs=NULL,
quiet=FALSE,
ssquiet=quiet,
... )
{
# preparation ------------------------------------------------------------------
StartTime <- Sys.time()
if(missing(datname)) datname <- deparse(substitute(dat))
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
# checks:
if( ! is.numeric(dat) ) stop("dat must be numeric.")
if(!is.vector(dat) & !quiet) on.exit(message("Note in distLfit: dat was not a vector."), add=TRUE)
if(length(gofProp)>1 | any(gofProp<0) | any(gofProp>1) )
  stop("gofProp must be a single value between 0 and 1.")
# remove NAs, convert to vector:
dat_full <- dat
dat <- as.numeric( dat[!is.na(dat)]  )
# truncate (fit values only to upper part of values):
dat <- dat[dat>=threshold]
#
# possible distributions: ------------------------------------------------------
dn <- lmomco::dist.list()
names(dn) <- dn
# Selection:
if( ! is.null(selection) )
  {
  if(!is.character(selection)) stop("Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector.")
  seldn <- !selection %in% dn
  if(any(seldn))
   {
   curseldn <- selection[seldn]
   if(!quiet) on.exit(message("Note in distLfit: selection (", toString(curseldn),
   ") not available in lmomco::dist.list(), thus removed."), add=TRUE)
   selection <- selection[!seldn]
   }
  dn <- dn[selection]
  names(dn) <- NULL
  }
else
# remove some to save time and errors, see ?dist.list # gld, gov and tri added
if(speed) dn <- dn[ ! dn %in%
   c("aep4","cau","emu","gep","gld","gov","kmu","kur","lmrq","sla","st3","texp","tri")]
# Check remaining sample size
if(length(dat) < 5) {if(!ssquiet)on.exit(message("Note in distLfit: sample size (",
                         length(dat), ") is too small to fit parameters (<5)."), add=TRUE)
  error_out <- as.list(dn) # this is very useful for distLquantile
  names(error_out) <- dn  # since it keeps the columns if a selection is given
  error_gof <- matrix(NA, nrow=length(dn), ncol=6)
  colnames(error_gof) <- c("RMSE", "R2", paste0("weight",1:3), "weightc")
  rownames(error_gof) <- dn
  return(list(dat=dat, datname=datname, gofProp=gofProp, parameter=error_out,
        gof=error_gof, error="dat size too small.",
        truncate=truncate, threshold=threshold, dat_full=dat_full))}
#
# Fit distribution parameters --------------------------------------------------
# L-Moments of sample  # package lmomco
mom <- lmomco::lmoms(dat, nmom=5)
# estimate parameters for each distribution:    # this takes time!
if(progbars) message("Parameter estimation from linear moments:")
parameter <- lapply(dn, function(d) lmomco::lmom2par(mom, type=d) )
# error catching:
if( length(parameter) != length(dn))
  {
  if(!quiet) on.exit(message("Note in distLfit: Some distributions could not be fitted. Possibly cau."), add=TRUE)
  names(parameter) <- sapply(parameter, "[[", "type")
  }
else names(parameter) <- dn
#
# Goodness of Fit, output list, plot -------------------------------------------
output <- distLgof(list(dat=dat, datname=datname, gofProp=gofProp, parameter=parameter,
                        truncate=truncate, threshold=threshold, dat_full=dat_full),
                   weightc=weightc, plot=FALSE, progbars=progbars, ks=ks, quiet=quiet)
# compare GOF
if(gofComp)
  {
  distLgofPlot(output, quiet=quiet)
  plot <- FALSE
  }
if(plot) output <- distLplot(dlf=output, cdf=cdf, legargs=legargs, histargs=histargs, ... )
if(!plot) output$coldist <- berryFunctions::rainbow2(if(is.null(selection)) 5 else length(selection))

if(time & !quiet) on.exit(message("distLfit execution took ", 
      signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds."), add=TRUE)
return(invisible(output))
} # end of function
