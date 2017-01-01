#' Fit distributions via L-moments
#' 
#' Fit several distributions via L-moments with \code{lmomco::\link[lmomco]{lmom2par}}
#' and compute goodness of fit measures.
#' 
#' @return invisible dlf object, see \code{\link{printL}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014, July 2015, Dec 2016
#' @seealso \code{\link{plotLfit}}, \code{\link{distLweights}}, \code{\link{plotLweights}}, 
#'     \code{extRemes::\link[extRemes]{fevd}}, \code{MASS::\link[MASS]{fitdistr}}.\cr
#'     More complex estimates of quality of fits:
#'     Fard, M.N.P. and Holmquist, B. (2013, Chilean Journal of Statistics): 
#'     Powerful goodness-of-fit tests for the extreme value distribution.
#'     http://chjs.mat.utfsm.cl/volumes/04/01/Fard_Holmquist(2013).pdf
#'          
#' @keywords hplot dplot distribution univar
#' @export
#' @importFrom lmomco dist.list lmoms lmom2par plmomco
#' @importFrom berryFunctions rainbow2 tryStack rmse rsquare
#' @importFrom utils getFromNamespace
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
#' # arguments that can be passed to plotting function:
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
#' y <- 10^rnorm(100, mean=2, sd=0.3) # if you use 1e4, distLfit will be much slower
#' hist(y, breaks=20)
#' berryFunctions::logHist(y, col=8)
#' dlf <- distLfit(log10(y))
#' plotLfit(dlf, breaks=50)
#' plotLfit(dlf, breaks=50, log=TRUE)
#' 
#' 
#' # Goodness of fit (see distLweights):
#' # measured by RMSE of cumulated distribution function and ?ecdf
#' 
#' 
#' 
#' # Fit all available distributions (30):
#' \dontrun{# this takes a while...
#' d_all <- distLfit(annMax, speed=FALSE) # 35 sec
#' printL(d_all)
#' plotLfit(d_all, nbest=22, distcols=grey(1:22/29), xlim=c(20,140))
#' plotLfit(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
#' plotLweights(dlf_all)
#' d_all$gof
#' }
#' 
#' @param dat       Vector with values
#' @param datname   Character string for main, xlab etc. 
#'                  DEFAULT: \code{deparse(substitute(dat))}
#' @param selection Selection of distributions. Character vector with types 
#'                  as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL
#' @param speed     If TRUE, several distributions are omitted, for the reasons 
#'                  shown in \code{lmomco::\link[lmomco]{dist.list}()}. DEFAULT: TRUE
#' @param ks        Include ks.test results and CDF R^2 in \code{dlf$gof}?
#'                  Computing is much faster when FALSE. DEFAULT: FALSE
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
#' @param \dots     Further arguments passed to \code{\link{distLweights}} 
#'                  like weightc, order=FALSE
#' 
distLfit <- function(
dat,
datname=deparse(substitute(dat)),
selection=NULL,
speed=TRUE,
ks=FALSE,
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
if(progbars) {  lapply <- pbapply::pblapply   ;   sapply <- pbapply::pbsapply  }
# checks:
if( ! is.numeric(dat) ) stop("dat must be numeric.")
if(!is.vector(dat) & !quiet) message("Note in distLfit: dat was not a vector.")
# remove NAs, convert to vector:
dat_full <- dat
dat <- as.numeric( dat[!is.na(dat)]  )
# truncate (fit values only to upper part of values):
dat <- dat[dat>=threshold] 
# GPD fits in q_gpd all use x>t, not x>=t, but if t=0, I want to use _all_ data
dat <- sort(dat, decreasing=TRUE)
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
#
# (error) output ---------------------------------------------------------------
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
if(lmomco::are.lmom.valid(mom))
{
  # estimate parameters for each distribution:    # this takes time!
  if(progbars) message("Parameter estimation from L-moments:")
  dlf$parameter <- lapply(dn, function(d) tryStack(lmomco::lmom2par(mom, type=d), silent=TRUE) )
} else 
{
  dlf$parameter <- rep(NA, length(dn))
  if(!quiet) message("Note in distLfit: L-moments are not valid. No distributions are fitted.")
  dlf$error <- c(error="invalid lmomco::lmoms", mom)
}
names(dlf$parameter) <- dn
#
# Error check ------------------------------------------------------------------
exclude <- sapply(dlf$parameter, function(x) 
  {
  if(is.null(x)) return(TRUE)
  if(inherits(x, "try-error")) return(TRUE)
  cumuprob <- try(lmomco::plmomco(mean(dlf$dat),x), silent=TRUE)
  if(is.null(cumuprob)||inherits(cumuprob,"try-error")) return(TRUE) 
  any(is.na(x$para))
  })
if(any(exclude))
  {
  dnexclude <- dn[exclude]
  if(!quiet) message("Note in distLfit: The following distributions were ",
                     "excluded since no parameters were estimated (",
                     sum(exclude),"/",length(dn),"):\n", toString(dnexclude),
                     if(sum(!exclude)<2) "\nGOF cannot be compared")
  dn <- dn[!exclude]
  # dlf$parameter <- dlf$parameter[!exclude] # not sure whether this is always good...
}
# Goodness of Fit --------------------------------------------------------------
# CDFS for RMSE (and R2): (dat must be sorted at this point in time!)
if(progbars) message("Calculating CDFs:")
tcdfs <- tryStack(lapply(dn, function(d) lmomco::plmomco(dat,dlf$parameter[[d]])))
names(tcdfs) <- dn # Theoretical CumulatedDensityFunctions
if(!quiet)
  {
  nNA <- base::sapply(tcdfs, function(x) sum(is.na(x)))
  if(any(nNA>0)) 
    message("Note in distLfit: there are NAs in CDF (distribution support region ",
            "probably does not span the whole data range): ", 
            toString(paste0(dn[nNA>0], " (", nNA[nNA>0], ")" )),
            " of ", length(tcdfs[[1]]), " values.")
  }
ecdfs <- ecdf(dlf$dat)(dat) # Empirical CDF
# RMSE:
if(progbars) message("Calculating RMSE:")
RMSE <- sapply(dn, function(d) berryFunctions::rmse(tcdfs[[d]], ecdfs, quiet=TRUE))
# add nonfitted distributions:
if(any(exclude)) 
  {
  RMSEexcl <- rep(NA, sum(exclude))
  names(RMSEexcl) <- dnexclude
  RMSE <- c(if(length(dn>0))RMSE,RMSEexcl)
  ###R2   <- c(  R2,RMSEexcl)
  }
# Weights for weighted averages:
dlf$gof <- distLweights(RMSE, ...)
# ks and R^2 values:
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) message("Performing ks.test:")
  ## library("lmomco") # flagged by R CMD check
  for(d in dn) assign(paste0("cdf",d), getFromNamespace(paste0("cdf",d), "lmomco"))
  ksA <- lapply(dn, function(d) ks.test(dlf$dat, paste0("cdf",d), dlf$parameter[[d]]) )
  ksP <- base::sapply(ksA, function(x) x$p.value   ) 
  ksD <- base::sapply(ksA, function(x) x$statistic )
  # R2
  if(progbars) message("Calculating R2:")
  R2 <- sapply(dn, function(d) berryFunctions::rsquare(tcdfs[[d]], ecdfs, quiet=TRUE))
  # add to output data.frame:
  dlf$gof$ksP[dn] <- ksP 
  dlf$gof$ksD[dn] <- ksD 
  dlf$gof$R2 [dn] <- R2 
  }
#
# time message + output --------------------------------------------------------
dlf$distnames <- rownames(dlf$gof)[order(dlf$gof$RMSE)] # in order even if order=FALSE
if(time & !quiet) message("distLfit execution took ", 
                  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(dlf))
} # end of function
