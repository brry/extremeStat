#' Fit distributions via L-moments
#' 
#' Fit several distributions via L-moments with \code{lmomco::\link[lmomco]{lmom2par}}.
#' Plot histogram and distribution densities \emph{or} ecdf with cumulated probability.
#' Compute goodness of fit values with \code{\link{distLgof}}.
#' 
#' @return invisible dlf object, see \code{\link{printL}}.
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
#' plotLfit(dlf)
#' 
#' # arguments that can be passed:
#' plotLfit(dlf, lty=2, col=3, legargs=list(lwd=3), main="booh!")
#' set.seed(42)
#' dlf_b <- distLfit(rbeta(100, 5, 2))
#' plotLfit(dlf_b, nbest=10, legargs=c(x="left"))
#' plotLfit(dlf_b, selection=c("gpa", "glo", "gev", "wak"))
#' plotLfit(dlf_b, selection=c("gpa", "glo", "gev", "wak"), order=TRUE)
#' plotLfit(dlf_b, distcols=c("orange",3:6), lty=1:3) # lty is recycled
#' plotLfit(dlf_b, cdf=TRUE)
#' plotLfit(dlf_b, cdf=TRUE, histargs=list(do.points=FALSE), sel="nor")
#' 
#' 
#' # logarithmic axes:
#' set.seed(1)
#' y <- 10^rnorm(100, mean=2, sd=0.3) # if you use 1e4, distLgof will be much slower
#' hist(y, breaks=20)
#' berryFunctions::logHist(y, col=8)
#' dlf <- distLfit(log10(y))
#' plotLfit(dlf, breaks=50)
#' plotLfit(dlf, breaks=50, log=TRUE)
#' 
#' \dontrun{
#' # this takes a while, as it tries to fit all 30 distributions:
#' d_all <- distLfit(annMax, speed=FALSE) # 35 sec
#' printL(d_all)
#' plotLfit(d_all, nbest=22, distcols=grey(1:22/29), xlim=c(20,140))
#' plotLfit(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
#' d_all$gof
#' }
#' 
#' @param dat       Vector with values
#' @param datname   Character string for main, xlab etc. 
#'                  DEFAULT: \code{deparse(substitute(dat))}
#' @param speed     If TRUE, several distributions are omitted, for the reasons 
#'                  shown in \code{lmomco::\link[lmomco]{dist.list}()}. DEFAULT: TRUE
#' @param selection Selection of distributions. Character vector with types 
#'                  as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL
#' @param truncate  Number between 0 and 1. POT Censored \code{\link{distLquantile}}: 
#'                  fit to highest values only (truncate lower proportion of x). 
#'                  Probabilities are adjusted accordingly. DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, 
#'                  set this only via truncate, see Details of \code{\link{q_gpd}}. 
#'                  DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param progbars  Show progress bars for each loop? DEFAULT: TRUE if n > 200
#' @param time      \code{\link{message}} execution time? DEFAULT: TRUE
#' @param quiet     Suppress notes? DEFAULT: FALSE
#' @param ssquiet   Suppress sample size notes? DEFAULT: quiet
#' @param \dots     Further arguments passed to \code{\link{distLgof}} 
#'                  like cweights, ks=TRUE, order=FALSE
#' 
distLfit <- function(
dat,
datname=deparse(substitute(dat)),
speed=TRUE,
selection=NULL,
truncate=0,
threshold=berryFunctions::quantileMean(dat, truncate),
progbars=length(dat)>200,
time=TRUE,
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
# possible / selected distributions --------------------------------------------
dn <- lmomco::dist.list()
names(dn) <- dn
# Selection:
if( ! is.null(selection) )
  {
  if(!is.character(selection)) stop("Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector.")
  seldn <- selection %in% dn
  if(any(!seldn))
   {
   if(!quiet) message("Note in distLfit: selection (", toString(selection[!seldn]),
                      ") not available in lmomco::dist.list(), thus removed.")
   selection <- selection[seldn]
   }
  dn <- dn[selection]
  }
else
# remove some to save time and errors, see ?dist.list # gld, gov and tri added
if(speed) dn <- dn[ ! dn %in%
   c("aep4","cau","emu","gep","gld","gov","kmu","kur","lmrq","sla","st3","texp","tri")]

# prepare output:
dlf <- list(dat=dat, dat_full=dat_full, datname=datname, 
            distnames=dn, distcols=berryFunctions::rainbow2(length(dn)), 
            distselector="distLfit", 
            truncate=truncate, threshold=threshold)

# Check remaining sample size:
if(length(dat) < 5) 
  {
  if(!ssquiet) message("Note in distLfit: sample size (", length(dat), 
                       ") is too small to fit parameters (<5).")
  # error output:
  dlf$parameter <- as.list(dn)
  dn2 <- rep(NA,length(dn)) ; names(dn2) <- dn
  dlf$gof <- suppressWarnings(distLweights(dn2))
  dlf$error <- paste0("dat size too small (",length(dat),")")
  return(invisible(dlf))
  }
#
# Fit distribution parameters --------------------------------------------------
# L-Moments of sample  # package lmomco
mom <- lmomco::lmoms(dat, nmom=5)
# estimate parameters for each distribution:    # this takes time!
if(progbars) message("Parameter estimation from L-moments:")
dlf$parameter <- lapply(dn, function(d) tryStack(lmomco::lmom2par(mom, type=d), silent=TRUE) )
names(dlf$parameter) <- dn
#
# Goodness of Fit --------------------------------------------------------------

dlf <- distLgof(dlf, progbars=progbars, quiet=quiet, ...)
dlf$distnames <- rownames(dlf$gof)[order(dlf$gof$RMSE)]

if(time & !quiet) message("distLfit execution took ", 
                  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(dlf))
} # end of function
