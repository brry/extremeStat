#' Quality of distribution fits
#' 
#' Calculate goodness of fit for several distributions, plot rank comparison.
#' 
#' @return List as explained in \code{\link{extremeStat}}. The added element is gof,\cr
#' a data.frame with the root mean square error (RMSE) and R squared (R2),\cr
#' if ks=TRUE, the p and D values from a simple ks.test,\cr
#' as well as weights by three different approaches for each distribution function.\cr
#' The weights are inverse to RMSE, weight1 for all dists, 
#' weight2 places zero weight on the worst function, weight3 on the worst half of functions.
#' @note If you get a \code{note in distLgof: NAs removed in CDF ...}, this
#'       probably means that the support of some of the fitted distributions do not
#'       span the whole data range. Instead the outside support regions get NAs that
#'       are then detected by rmse and rsquare. I plan to fix this with WHA's new supdist.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014 + July 2015
#' @seealso \code{\link{distLweights}}, \code{\link{plotLgof}}, \code{\link{distLfit}} 
#'     More complex estimates of quality of fits:
#'     Fard, M.N.P. and Holmquist, B. (2013, Chilean Journal of Statistics): 
#'     Powerful goodness-of-fit tests for the extreme value distribution.
#'     http://chjs.mat.utfsm.cl/volumes/04/01/Fard_Holmquist(2013).pdf
#' @keywords univar hplot distribution
#' @export
#' @importFrom lmomco plmomco
#' @importFrom berryFunctions rmse rsquare
#' @importFrom utils getFromNamespace
#' 
#' @examples
#' 
#' library(lmomco)
#' data(annMax)
#' 
#' # Goodness of Fit is measured by RMSE of cumulated distribution function and ?ecdf
#' dlf <- distLfit(annMax, cdf=TRUE, nbest=17)
#' plotLfit(dlf, cdf=TRUE, sel=c("wak", "revgum"))
#' dlf$gof
#' x <- sort(annMax)
#' segments(x0=x, y0=plmomco(x, dlf$parameter$revgum), y1=ecdf(annMax)(x), col=2)
#' segments(x0=x, y0=plmomco(x, dlf$parameter$wak), y1=ecdf(annMax)(x), col=4, lwd=2)
#' plot(ecdf(annMax), add=TRUE)
#' # RMSE: root of average of ( errors squared )  ,   errors = line distances
#' 
#' # weights by three different weighting schemes
#' plotLgof(dlf)
#' distLgof(dlf, ks=TRUE, plot=FALSE)$gof
#' 
#' # Kolmogorov-Smirnov Tests return slightly different values:
#' ks.test(annMax, "pnorm", mean(annMax), sd(annMax) )$p.value
#' ks.test(annMax, "cdfnor", parnor(lmoms(annMax)))$p.value
#' 
#' 
#' # GOF: how well do the distributions fit the original data?
#' pairs(dlf$gof[,1:3]) # measures of goodness of fit are correlated quite well here.
#' 
#' 
#' \dontrun{ ## to save CRAN check computing time
#' 
#' plotLfit(dlf, cdf=TRUE, sel=c("pe3", "rice", "revgum"), order=T)
#' x <- sort(annMax, decreasing=TRUE)[  1:(0.5*length(annMax))  ]
#' tcdfs <- plmomco(x,dlf$parameter[["revgum"]])
#' ecdfs <- ecdf(annMax)(x) # Empirical CDF
#' plot(x, tcdfs, type="o", col=2)
#' points(x, ecdfs)
#' linReg(tcdfs, ecdfs, type="o")
#' abline(a=0,b=1, lty=3)
#' 
#' } # end dontrun
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the elements
#'            \code{dat, datname, parameter}
#' @param order Sort output$gof by RMSE? If FALSE, the order of appearance in 
#'              selection (or dlf$parameter) is kept. DEFAULT: TRUE
#' @param plot Call \code{\link{plotLgof}}? DEFAULT: TRUE
#' @param progbars Show progress bars for each loop? DEFAULT: TRUE if n > 200
#' @param ks Include ks.test results in \code{dlf$gof}?
#'            Computing is much faster when FALSE. DEFAULT: TRUE
#' @param weightc Custom weights, see \code{\link{distLweights}}. DEFAULT: NA
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param \dots Further arguments passed to \code{\link{plotLgof}}
#' 
distLgof <- function(
dlf,
order=TRUE,
plot=TRUE,
progbars=length(dlf$dat)>200,
ks=TRUE,
weightc=NA,
quiet=FALSE,
...
)
{
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
dn <- names(dlf$parameter)
# Error check:
exclude <- sapply(dlf$parameter, function(x) 
  {
  if(is.null(x)) return(TRUE)
  # if("ifail" %in% names(x)) if(x$ifail != 0) return(TRUE) ## restriction too tight
  ### if(x$type=="gld") return(TRUE) # lmomco 2.2.2 cdfgld bug
  if(inherits(x, "try-error")) return(TRUE) # lmomco <=2.2.4: parkap TAU4 NA error
  cumuprob <- try(lmomco::plmomco(mean(dlf$dat),x), silent=TRUE)
  if(is.null(cumuprob)||inherits(cumuprob,"try-error")) return(TRUE)  # kappa errors
  any(is.na(x$para))
  })
if(any(exclude))
  {
  dnexclude <- dn[exclude]
  if(!quiet) message("Note in distLgof: The following distributions were ",
                     "excluded since no parameters were estimated:\n",
                      toString(dnexclude))
  dn <- dn[!exclude]
  # dlf$parameter <- dlf$parameter[!exclude] # not sure whether this is always good...
}
if(length(dn)<1&!quiet) message("Note in distLgof: No fitted distributions",
                               " in dlf, thus GOF can't be compared.") else
if(length(dn)<2&!quiet) message("Note in distLgof: Only ", toString(dn),
                                " was fitted, thus GOF can't be compared.")
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) message("Performing ks.test:")
  ## library("lmomco") # flagged by R CMD check
  for(d in dn) assign(paste0("cdf",d), getFromNamespace(paste0("cdf",d), "lmomco"))
  ksA <- lapply(dn, function(d) ks.test(dlf$dat, paste0("cdf",d), dlf$parameter[[d]]) )
  ksP <- sapply(ksA, function(x) x$p.value   )
  ksD <- sapply(ksA, function(x) x$statistic )
  names(ksD) <- dn
  }
# CDFS for R2 and RMSE:
dat2 <- sort(dlf$dat, decreasing=TRUE)
if(progbars) message("Calculating CDFs:")
tcdfs <- lapply(dn, function(d) lmomco::plmomco(dat2,dlf$parameter[[d]]))
names(tcdfs) <- dn # Theoretical CumulatedDensityFunctions
ecdfs <- ecdf(dlf$dat)(dat2) # Empirical CDF
# Root Mean Square Error, R squared:
if(progbars) sapply <- pbapply::pbsapply
if(progbars) message("Calculating RMSE:")
RMSE <- sapply(dn, function(d)    berryFunctions::rmse(tcdfs[[d]], ecdfs, quiet=TRUE))
if(progbars) message("Calculating R2:")
R2   <- sapply(dn, function(d) berryFunctions::rsquare(tcdfs[[d]], ecdfs, quiet=TRUE))
if(!quiet)
  {
  nNA <- base::sapply(tcdfs, function(x) sum(is.na(x)))
  if(any(nNA>0)) 
    {
    dNA <- paste(paste0(dn[nNA>0], " (", nNA[nNA>0], ")"), collapse=", ")
    message("Note in distLgof: NAs removed in CDF (limited support region?): ", 
                    dNA, " of ", length(tcdfs[[1]]), " values.")
    }
  }
# if all distributions are excluded:
if(length(RMSE)==0 && FALSE) # should be obsolete with next part
  {
  RMSE <- rep(NA, length(dlf$parameter)) 
  R2   <- rep(NA, length(dlf$parameter)) 
  names(RMSE) <- names(dlf$parameter)
  }
# add nonfitted distributions:
if(any(exclude)) 
  {
  RMSEexcl <- rep(NA, sum(exclude))
  names(RMSEexcl) <- dnexclude
  RMSE <- c(RMSE,RMSEexcl)
  R2   <- c(  R2,RMSEexcl)
  }
# Weights for weighted averages:
gof <- data.frame(RMSE, R2)
if(ks) {gof$ksP <- ksP; gof$ksD <- ksD}
if(order) gof <- gof[ order(RMSE), ]
gof <- cbind(gof, distLweights(RMSE, order=order, weightc=weightc))
# output:
dlf$gof <- gof
if(plot) plotLgof(dlf, quiet=quiet, ...)
dlf
} # end of function

