#' Fit distributions via linear moments
#' 
#' Fit several distributions via linear moments, plot histogram and
#' distribution densities \emph{or} ecdf with cumulated probability.
#' Also returns goodness of fit values.
#' This is the main fitting function calling   distLgof
#'
#' @details Fits parameters via \code{\link[lmomco]{lmom2par}} in the package \code{lmomco}
#' @return List as explained in \code{\link{extremeStat}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014 + July 2015
#' @seealso \code{\link{distLgof}}, \code{\link{plotLfit}}.
#'          \code{\link[extRemes]{fevd}} in the package \code{extRemes},
#'          \code{\link[MASS]{fitdistr}} in the package \code{MASS}.
#' @keywords hplot dplot distribution
#' @export
#' @importFrom lmomco dist.list lmoms lmom2par
#' @importFrom berryFunctions rainbow2 tryStack
#' 
#' @examples
#' 
#' data(annMax)
#' # basic usage on real data (annual discharge maxima in Austria)
#' dlf <- distLfit(annMax)
#' str(dlf, max.lev=2)
#' printL(dlf)
#' 
#' # arguments that can be passed:
#' distLfit(annMax, lty=2, col=3, legargs=list(lwd=3), main="booh!")
#' set.seed(42)
#' dlf_b <- distLfit(rbeta(100, 5, 2), nbest=10, legargs=c(x="left"))
#' plotLfit(dlf_b, selection=c("gpa", "glo", "gev", "wak"))
#' plotLfit(dlf_b, selection=c("gpa", "glo", "gev", "wak"), order=TRUE)
#' plotLfit(dlf_b, coldist=c("orange",3:6), lty=1:3) # lty is recycled
#' plotLfit(dlf_b, cdf=TRUE)
#' plotLfit(dlf_b, cdf=TRUE, histargs=list(do.points=FALSE), sel="nor")
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
#' plotLfit(dlf, breaks=50, log=TRUE)
#' 
#' \dontrun{
#' # this takes a while, as it tries to fit all 30 distributions:
#' d_all <- distLfit(annMax, speed=FALSE, plot=FALSE) # 35 sec
#' printL(d_all)
#' plotLfit(d_all, nbest=22, coldist=grey(1:22/29), xlim=c(20,140))
#' plotLfit(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
#' d_all$gof
#' }
#' 
#' @param dat       Vector with values
#' @param datname   Character string for main, xlab etc. 
#'                  DEFAULT: \code{deparse(substitute(dat))}
#' @param speed     If TRUE, several distributions are omitted, for the reasons 
#'                  shown in \code{lmomco::\link[lmomco]{dist.list}()}. DEFAULT: TRUE
#' @param ks        Include ks.test results in \code{dlf$gof}? Computing is much 
#'                  faster when FALSE. DEFAULT: FALSE
#' @param selection Selection of distributions. Character vector with types 
#'                  as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL
#' @param order     Should gof output be ordered by fit? 
#'                  Strongly recommended with the default color palette. DEFAULT: TRUE
#' @param weightc   Named custom weights for each distribution, 
#'                  see \code{\link{distLgof}}. DEFAULT: NA
#' @param truncate  Number between 0 and 1. POT Censored \code{\link{distLquantile}}: 
#'                  fit to highest values only (truncate lower proportion of x). 
#'                  Probabilities are adjusted accordingly. DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, 
#'                  set this only via truncate, see Details of \code{\link{q_gpd}}. 
#'                  DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param progbars  Show progress bars for each loop? DEFAULT: TRUE if n > 200
#' @param time      \code{\link{message}} execution time? DEFAULT: TRUE
#' @param plot      Should a histogram with densities be plotted? DEFAULT: TRUE
#' @param cdf       If TRUE, plot cumulated DF instead of probability density. 
#'                  DEFAULT: FALSE
#' @param legargs   List of arguments passed to \code{\link{legend}} 
#'                  except for legend and col. DEFAULT: NULL
#' @param histargs  List of arguments passed to \code{\link{hist}} 
#'                  except for x, breaks, col, xlim, freq. DEFAULT: NULL
#' @param quiet     Suppress notes? DEFAULT: FALSE
#' @param ssquiet   Suppress sample size notes? DEFAULT: quiet
#' @param \dots     Further arguments passed to \code{\link{plotLfit}} 
#'                  if they are accepted there, else passed to 
#'                  \code{\link{lines}}, like lty, type, pch, ...
#' 
distLfit <- function(
dat,
datname=deparse(substitute(dat)),
speed=TRUE,
ks=FALSE,
selection=NULL,
order=TRUE,
weightc=NA,
truncate=0,
threshold=berryFunctions::quantileMean(dat, truncate),
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
datname <- datname # evaluate promise before dat is changed
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
# checks:
if( ! is.numeric(dat) ) stop("dat must be numeric.")
if(!is.vector(dat) & !quiet) message("Note in distLfit: dat was not a vector.")
# remove NAs, convert to vector:
dat_full <- dat
dat <- as.numeric( dat[!is.na(dat)]  )
# truncate (fit values only to upper part of values):
dat <- dat[dat>=threshold] # GPD fits in q_gpd all use x>t, not x>=t
#                          # but if t=0, I want to use _all_ data
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
   if(!quiet) message("Note in distLfit: selection (", toString(selection[seldn]),
                      ") not available in lmomco::dist.list(), thus removed.")
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
if(length(dat) < 5) {if(!ssquiet) message("Note in distLfit: sample size (",
                                  length(dat), ") is too small to fit parameters (<5).")
  error_out <- as.list(dn) # this is very useful for distLquantile
  names(error_out) <- dn  # since it keeps the rows if a selection is given
  error_gof <- matrix(NA, nrow=length(dn), ncol=6)
  colnames(error_gof) <- c("RMSE", "R2", paste0("weight",1:3), "weightc")
  rownames(error_gof) <- dn
  return(list(dat=dat, datname=datname, parameter=error_out,
        gof=error_gof, error="dat size too small.",
        truncate=truncate, threshold=threshold, dat_full=dat_full))}
#
# Fit distribution parameters --------------------------------------------------
# L-Moments of sample  # package lmomco
mom <- lmomco::lmoms(dat, nmom=5)
# estimate parameters for each distribution:    # this takes time!
if(progbars) message("Parameter estimation from linear moments:")
parameter <- lapply(dn, function(d) tryStack(lmomco::lmom2par(mom, type=d), silent=TRUE) )
# wrapped in try since july 2016 because parkap breaks if TAU4=NA  (lmomco 2.2.4)
# error catching:
if( length(parameter) != length(dn))
  {
  if(!quiet) message("Note in distLfit: Some distributions could not be fitted. Possibly cau.")
  names(parameter) <- sapply(parameter, "[[", "type")
  }
else names(parameter) <- dn
#
# Goodness of Fit, output list, plot -------------------------------------------
output <- distLgof(list(dat=dat, datname=datname, parameter=parameter,
                        truncate=truncate, threshold=threshold, dat_full=dat_full),
     weightc=weightc, plot=FALSE, progbars=progbars, ks=ks, quiet=quiet, order=order)

if(plot) output <- plotLfit(dlf=output, cdf=cdf, legargs=legargs, histargs=histargs, ... )
if(!plot) output$coldist <- berryFunctions::rainbow2(if(is.null(selection)) 5 else length(selection))

if(time & !quiet) message("distLfit execution took ", 
                  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(output))
} # end of function
