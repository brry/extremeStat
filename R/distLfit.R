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
#' plotLfit(dlf, lty=2, col=3, nbest=17, legargs=list(lwd=3), main="booh!")
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
#' y <- 10^rnorm(300, mean=2, sd=0.3) # if you use 1e4, distLfit will be much slower
#' hist(y, breaks=20)
#' berryFunctions::logHist(y, col=8)
#' dlf <- distLfit(log10(y))
#' plotLfit(dlf, breaks=50)
#' plotLfit(dlf, breaks=50, log=TRUE)
#' 
#' 
#' # Goodness of fit: how well do the distributions fit the original data?
#' # measured by RMSE of cumulated distribution function and ?ecdf
#' # RMSE: root of average of ( errors squared )  ,   errors = line distances
#' dlf <- distLfit(annMax, ks=TRUE)
#' plotLfit(dlf, cdf=TRUE, sel=c("wak", "revgum"))
#' x <- sort(annMax)
#' segments(x0=x, y0=lmomco::plmomco(x, dlf$parameter$revgum), y1=ecdf(annMax)(x), col=2)
#' segments(x0=x, y0=lmomco::plmomco(x, dlf$parameter$wak), y1=ecdf(annMax)(x), col=4, lwd=2)
#' # weights by three different weighting schemes, see distLweights:
#' plotLweights(dlf)
#' plotLfit(distLfit(annMax              ), cdf=TRUE, nbest=17)$gof
#' plotLfit(distLfit(annMax, truncate=0.7), cdf=TRUE, nbest=17)$gof
#' pairs(dlf$gof[,-(2:5)]) # measures of goodness of fit are correlated quite well here.
#' dlf$gof
#'
#' # Kolmogorov-Smirnov Tests for normal distribution return slightly different values:
#' library(lmomco)
#' ks.test(annMax, "pnorm", mean(annMax), sd(annMax) )$p.value
#' ks.test(annMax, "cdfnor", parnor(lmoms(annMax)))$p.value
#' 
#' 
#' # Fit all available distributions (30):
#' \dontrun{# this takes a while...
#' d_all <- distLfit(annMax, speed=FALSE, progbars=TRUE) # 20 sec
#' printL(d_all)
#' plotLfit(d_all, nbest=30, distcols=grey(1:22/29), xlim=c(20,140))
#' plotLfit(d_all, nbest=30, ylim=c(0,0.04), xlim=c(20,140))
#' plotLweights(d_all)
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
if(progbars) lapply <- pbapply::pblapply
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
# remove some to save time and errors, see ?dist.list # gld, gov and tri added
if(speed) dn <- dn[ ! dn %in%
   c("aep4","cau","emu","gep","gld","gov","kmu","kur","lmrq","sla","st3","texp","tri")]
#
# Selection:
if( ! is.null(selection) )
  {
  if(!is.character(selection)) stop("Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector.")
  nondn <- selection[!selection %in% dn]
  if(length(nondn)>0 & !quiet) message("Note in distLfit: selection (", toString(nondn),
                                       ") not available in lmomco::dist.list().")
  dn <- selection ; names(dn) <- dn
  }
#
# output list ------------------------------------------------------------------
dlf <- list(parameter=as.list(replace(dn,TRUE,NA)),
            dat_full=dat_full, dat=dat, datname=datname, 
            distnames=dn, distcols=berryFunctions::rainbow2(length(dn)), 
            distselector="distLfit", distfailed="",
            truncate=truncate, threshold=threshold)

# Check remaining sample size:
if(length(dat) < 5) 
  {
  if(!ssquiet) message("Note in distLfit: sample size (", length(dat), 
                       ") is too small to fit parameters (<5).")
  # error output:
  dlf$gof <- distLweights(replace(dn,TRUE,NA), quiet=TRUE, ...)
  dlf$error <- paste0("dat size too small (",length(dat),")")
  dlf$distfailed <- dn
  return(invisible(dlf))
  }
#
# Fit distribution parameters --------------------------------------------------
# This is the actual work...
# L-Moments of sample  # package lmomco
mom <- lmomco::lmoms(dat, nmom=5)
if(lmomco::are.lmom.valid(mom))
{
  # estimate parameters for each distribution:    # this takes time!
  if(progbars) message("Parameter estimation from L-moments:")
  dlf$parameter <- lapply(dn, function(d) tryStack(lmomco::lmom2par(mom, type=d), silent=TRUE) )
} else 
{
  if(!quiet) message("Note in distLfit: L-moments are not valid. No distributions are fitted.")
  dlf$gof <- distLweights(replace(dn,TRUE,NA), quiet=TRUE, ...)
  dlf$error <- c(error="invalid lmomco::lmoms", mom)
  dlf$distfailed <- dn
}
names(dlf$parameter) <- dn
#
# Error check ------------------------------------------------------------------
failed <- sapply(dlf$parameter, function(x) 
  {
  if(is.null(x)) return(TRUE)
  if(all(is.na(x))) return(TRUE)
  if(inherits(x, "try-error")) return(TRUE)
  if(x$type=="kap") 
    {
    if(x$ifail!=0) return(TRUE)
    if(round(x$support[2],7)<=round(x$support[1],7) ) return(TRUE)
    quant <- lmomco::quakap(seq(truncate,1,len=200), para=x) 
    # plot(xxx<-seq(-100,50,len=900), lmomco::pdfkap(xxx, para=x), type="l")
    if(length(unique(quant))<20) return(TRUE) # xx5 in test-quantile
    }
  cumuprob <- suppressWarnings(try(lmomco::plmomco(mean(dlf$dat),x), silent=TRUE))
  if(is.null(cumuprob)||inherits(cumuprob,"try-error")) return(TRUE) 
  any(is.na(x$para))
  })
if(any(failed))
  {
  dlf$parameter[failed] <- NA # needed in kappa cases like xx3 in test-quantile
  dlf$distfailed <- dn[failed]
  if(!quiet) message("Note in distLfit: ", sum(failed),"/",length(dn),
                     " distributions could not be fitted: ", toString(dlf$distfailed),
                     if(sum(!failed)<2) "\nGOF cannot be compared")
  ###dn <- dn[!failed]
}
# Goodness of Fit --------------------------------------------------------------
# CDFS for RMSE (and R2): (dat must be sorted at this point)
if(progbars) message("Calculating CDFs:")
# Theoretical CumulatedDensityFunctions:
tcdfs <- suppressWarnings(
         lapply(dn, function(d) tryStack(lmomco::plmomco(dat,dlf$parameter[[d]]), silent=TRUE)))
# support region check:
if(!quiet)
  {
  nNA <- sapply(tcdfs, function(x) sum(is.na(x)))
  if(any(nNA>0)) 
    message("Note in distLfit: there are NAs in CDF (distribution support region ",
            "probably does not span the whole data range): ", 
            toString(paste0(dn[nNA>0], " (", nNA[nNA>0], ")" )),
            " of ", length(tcdfs[[1]]), " values.")
  }
ecdfs <- ecdf(dlf$dat)(dat) # Empirical CDF
# rescale for truncation
tcdfs[sapply(tcdfs, inherits, "try-error")] <- NA
# RMSE:
if(progbars) message("Calculating RMSE:")
RMSE <- suppressWarnings(
        lapply(dn, function(d)  tryStack(rmse(tcdfs[[d]], ecdfs, quiet=TRUE), 
                                         silent=TRUE)))
# change nonfitted distributions RMSE to NA:
RMSE <- sapply(RMSE, function(x) if(inherits(x, "try-error")) NA else x)
# distribution weights:
dlf$gof <- distLweights(RMSE, quiet=if(sum(!failed)<2) TRUE else quiet, ...)
# ks and R^2 values:
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) message("Performing ks.test:")
  ## library("lmomco") # flagged by R CMD check, not necessary if plmomco is imported
  ksA <- suppressWarnings(
         lapply(dn, function(d) tryStack(ks.test(dlf$dat, "plmomco", 
                                                 dlf$parameter[[d]]), silent=TRUE)))
  ksP <- sapply(ksA, function(x) x$p.value   ) 
  ksD <- sapply(ksA, function(x) x$statistic )
  # R2
  if(progbars) message("Calculating R2:")
  R2 <- suppressWarnings(
        lapply(dn, function(d)  tryStack(rsquare(tcdfs[[d]], ecdfs, quiet=TRUE), 
                                         silent=TRUE)))
  R2 <- sapply(R2, function(x) if(inherits(x, "try-error")) NA else x)
  # add to output data.frame:
  dlf$gof$ksP <- ksP[rownames(dlf$gof)] 
  dlf$gof$ksD <- ksD[paste0(rownames(dlf$gof),".D")]
  dlf$gof$R2  <-  R2[rownames(dlf$gof)] 
  }
#
# time message + output --------------------------------------------------------
dlf$distnames <- rownames(dlf$gof)
if(time & !quiet) message("distLfit execution took ", 
                  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.")
return(invisible(dlf))
} # end of function
