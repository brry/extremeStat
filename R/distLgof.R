#' Quality of distribution fits
#' 
#' Calculate goodness of fit for several distributions, plot rank comparison.
#' 
#' @return List as explained in \code{\link{extremeStat}}. The added element is gof,\cr
#' a data.frame the root mean square error (RMSE) and R squared (R2) for the top \code{gofProp} of \code{dat},\cr
#' if ks=TRUE, the p and D values from a simple ks.test,\cr
#' as well as weights by three different approaches for each distribution function.\cr
#' The weights are inverse to RMSE, weight1 for all dists, weight2 places zero weight on the worst function, weight3 on the worst half of functions.
#' @note If you get a \code{note in distLgof: NAs removed in CDF ...}, this
#'       probably means that the support of some of the fitted distributions do not
#'       span the whole data range. Instead the outside support regions get NAs that
#'       are then detected by rmse and rsquare. I plan to fix this with WHA's new supdist.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014 + July 2015
#' @seealso \code{\link{distLgofPlot}}, \code{\link{distLfit}}. More complex estimates of quality of fits:
#'          \url{http://chjs.soche.cl/papers/vol4n1_2013/ChJS-04-01-04.pdf}
#' @keywords univar hplot distribution
#' @export
#' @importFrom lmomco plmomco
#' @importFrom berryFunctions rmse rsquare
#' 
#' @examples
#' 
#' library(lmomco)
#' data(annMax)
#' 
#' # Goodness of Fit is measured by RMSE of cumulated distribution function and ?ecdf
#' dlf <- distLfit(annMax, cdf=TRUE, nbest=17)
#' distLplot(dlf, cdf=TRUE, sel=c("wak", "revgum"))
#' dlf$gof
#' x <- sort(annMax)
#' segments(x0=x, y0=plmomco(x, dlf$parameter$revgum), y1=ecdf(annMax)(x), col=2)
#' segments(x0=x, y0=plmomco(x, dlf$parameter$wak), y1=ecdf(annMax)(x), col=4, lwd=2)
#' plot(ecdf(annMax), add=TRUE)
#' # RMSE: root of average of ( errors squared )  ,   errors = line distances
#' 
#' # weights by three different weighting schemes
#' distLgofPlot(dlf, ranks=FALSE, weights=TRUE)
#' distLgof(dlf, ks=TRUE, plot=FALSE)$gof
#' 
#' # Kolmogorov-Smirnov Tests return slightly different values:
#' ks.test(annMax, "pnorm", mean(annMax), sd(annMax) )$p.value
#' ks.test(annMax, "cdfnor", parnor(lmoms(annMax)))$p.value
#' 
#' 
#' # GOF: how well do the distributions fit the original data?
#' pairs(dlf$gof[,1:3]) # measures of goodness of fit are correlated quite well here.
#' # In the next example, however, we see that it does matter which one is used.
#' 
#' # compare Ranks with different proportions used for GOF
#' par(mfrow=c(1,2))
#' distLgofPlot(dlf, weights=FALSE)
#' d <- distLgof(dlf, gofProp=0.8, plot=TRUE, weights=FALSE)
#' par(mfrow=c(1,1))
#' 
#' 
#' # effect of Proportion of values used to calculate RMSE
#' dlf100 <- distLfit(annMax, gofProp=1, nbest=19, breaks=10) # the default gofProp: 100%
#' distLgofPlot(dlf100)
#' 
#' dlf50 <- distLgof(dlf100, gofProp=0.5)
#' # so revgum, nor and rice do well on the upper half by R2, but bad by RMSE
#' distLplot(dlf50, breaks=10)
#' # The red dashed line shows the cut above which the data were used to get R2/rmse
#' 
#' distLplot(dlf=dlf50, cdf=TRUE)
#' distLplot(dlf=dlf50, selection=c("pe3","wei", "rice", "nor", "revgum"),
#'           xlim=c(60,120), ylim=c(0.5, 1), cdf=TRUE, col=1)
#' dlf50$gof
#' 
#' compranks <- function(d)
#' {
#' gofProp <- 0.5
#' x <- sort(annMax, decreasing=TRUE)[  1:(gofProp*length(annMax))  ]
#' tcdfs <- plmomco(x,dlf50$parameter[[d]])
#' ecdfs <- ecdf(annMax)(x) # Empirical CDF
#' # Root Mean Square Error, R squared:
#' berryFunctions::linReg(tcdfs, ecdfs, lwd=1, pch=16, main=d, digits=5, xlim=c(0.5, 1),
#'        ylim=c(0.5, 1), pos1="topleft")
#' abline(a=0,b=1, lty=3)
#' c(berryFunctions::rmse(tcdfs, ecdfs), berryFunctions::rsquare(tcdfs, ecdfs))
#' }
#' dn <- rownames(dlf50$gof)
#' 
#' op <- par(mfrow=c(5,4), mar=rep(1,4), xaxt="n", yaxt="n")
#' for(i in dn) compranks(i)
#' par(op)
#' # so revgum, nor and rice systematically deviate from ECDF.
#' # RMSE is better to sort by than R2
#' 
#' 
#' \dontrun{ ## to save CRAN check computing time
#' 
#' # custom weights
#' cw <- c("gpa"=7, "gev"=3, "wak"=6, "wei"=4, "kap"=3.5, "gum"=3, "ray"=2.1,
#'         "ln3"=2, "pe3"=2.5, "gno"=4, "gam"=5)
#' dlf <- distLfit(annMax, plot=FALSE, weightc=cw)
#' distLgofPlot(dlf, ranks=TRUE)
#' 
#' 
#' dev.new()
#' distLplot(dlf50, cdf=TRUE, sel=c("pe3", "rice", "revgum"), order=T)
#' x <- sort(annMax, decreasing=TRUE)[  1:(0.5*length(annMax))  ]
#' tcdfs <- plmomco(x,dlf50$parameter[["revgum"]])
#' ecdfs <- ecdf(annMax)(x) # Empirical CDF
#' plot(x, tcdfs, type="o", col=2)
#' points(x, ecdfs)
#' dev.new()
#' linReg(tcdfs, ecdfs, type="o")
#' abline(a=0,b=1, lty=3)
#' 
#' dev.new(record=TRUE)
#' for(i in 1:10/10) distLgof(dlf100, gofProp=i, weights=FALSE,
#'                            main=paste("upper", i*100, "% used for R2"))
#' # depending on which proportion of the data the GOF is calculated with, different
#' # distributions are selected to be the 5 best.
#' dev.off()
#' } # end dontrun
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the elements
#'            \code{dat, datname, parameter, gofProp}
#' @param gofProp Overrides value in list. Proportion (0:1) of highest values in
#'                \code{dat} to compute goodness of fit (dist / ecdf) with.
#'                 This enables to focus on the dist tail
#' @param plot Call \code{\link{distLgofPlot}}? DEFAULT: TRUE
#' @param progbars Show progress bars for each loop? DEFAULT: TRUE if n > 200
#' @param ks Include ks.test results in \code{dlf$gof}?
#'            Computing is much faster when FALSE. DEFAULT: TRUE
#' @param weightc Optional: a named vector with custom weights for each distribution.
#'                Are internally normalized to sum=1 after removing nonfitted dists.
#'                Names must match the parameter names from \code{\link{distLfit}}.
#'                DEFAULT: NA
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param \dots Further arguments passed to \code{\link{distLgofPlot}}
#' 
distLgof <- function(
dlf,
gofProp,
plot=TRUE,
progbars=length(dlf$dat)>200,
ks=TRUE,
weightc=NA,
quiet=FALSE,
...
)
{
if(any(!is.na(weightc))) if(is.null(names(weightc))) stop("weightc must have names.")
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
# gofProp overwrite and check:
if(!missing(gofProp)) dlf$gofProp <- gofProp
if(length(dlf$gofProp)>1 | any(dlf$gofProp<0) | any(dlf$gofProp>1) ) stop("gofProp must be a single value between 0 and 1.")
dn <- names(dlf$parameter)
# Error check:
exclude <- sapply(dlf$parameter, function(x) 
  {
  if(is.null(x)) return(TRUE)
  # if("ifail" %in% names(x)) if(x$ifail != 0) return(TRUE) ## restriction too tight
  ### if(x$type=="gld") return(TRUE) # lmomco 2.2.2 cdfgld bug
  if(inherits(x, "try-error")) return(TRUE) # lmomco <=2.2.4: parkap TAU4 NA error
  if(is.null(lmomco::plmomco(mean(dlf$dat),x))) return(TRUE)  # *)
  any(is.na(x$para))
  })                    #  *): CDF cannot be computed for kappa in Dresden example
if(any(exclude))
  {
  curdnexclude <- dn[exclude]
  if(!quiet) on.exit(message("Note in distLgof: The following distributions were excluded since no parameters were estimated:\n",
             toString(curdnexclude)), add=TRUE)
  dn <- dn[!exclude]
  # dlf$parameter <- dlf$parameter[!exclude] # not sure whether this is always good...
}
if(length(dn)<1&!quiet) on.exit(message("Note in distLgof: No fitted distributions",
                               " in dlf, thus GOF can't be compared."), add=TRUE) else
if(length(dn)<2&!quiet) on.exit(message("Note in distLgof: Only ", toString(dn),
                                " was fitted, thus GOF can't be compared."), add=TRUE)
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) message("Performing ks.test:")
  ksA <- lapply(dn, function(d) ks.test(dlf$dat, paste0("cdf",d), dlf$parameter[[d]]) )
  ksP <- sapply(ksA, function(x) x$p.value   )
  ksD <- sapply(ksA, function(x) x$statistic )
  names(ksD) <- dn
  }
# CDFS for R2 on upper gofProp of data:
dat2 <- sort(dlf$dat, decreasing=TRUE)[  1:(dlf$gofProp*length(dlf$dat))  ]
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
    on.exit(message("Note in distLgof: NAs removed in CDF (limited support region?): ", 
                    dNA, " of ", length(tcdfs[[1]]), " values."), add=TRUE)
    }
  }
# All into one data.frame:
gof <- data.frame(RMSE=RMSE, R2=R2)
if(all(dim(gof) == 0)) # dim = 0,0 if all distributions are excluded
{
gof <- data.frame(matrix(NA, ncol=5, nrow=length(dlf$parameter) ))
colnames(gof) <- c("RMSE","R2","weight1","weight2","weight3")
rownames(gof) <-  names(dlf$parameter)
} else
{
if(ks) {gof$ksP=ksP; gof$ksD=ksD}
# order by GOF:
gof <- gof[ order(gof$RMSE), ]  # sorting by gof$R2 does not work, see examples compranks
# Weights for weighted average of return values:
# Weight including all distributions
gof$weight1 <- gof[,"RMSE"] # the lower, the better, the more weight
gof$weight1 <- max(gof$weight1)-gof$weight1+min(gof$weight1)  # with min or mean added,
gof$weight1 <- gof$weight1/sum(gof$weight1)  # the worst fit is not completely excluded
# Exclude worst fit (needs 2 or more distributions in selection to work)
gof$weight2 <- gof[,"RMSE"]
gof$weight2 <- max(gof$weight2)-gof$weight2
gof$weight2 <- gof$weight2/sum(gof$weight2)
# use only best half (needs 4 or more dists - technically, 3 should be enough...)
gof$weight3 <- gof[,"RMSE"]
gof$weight3 <- max(gof$weight3)-gof$weight3
gof$weight3[(nrow(gof)/2):nrow(gof)] <- 0
gof$weight3 <- gof$weight3/sum(gof$weight3)
# custom weight
gof$weightc <- 0
if(any(!is.na(weightc)))
  {
  weightc <- weightc[names(weightc)%in%rownames(gof)]
  gof[names(weightc), "weightc"] <- weightc
  }
gof$weightc <- gof$weightc/sum(gof$weightc)
# add nonfitted distributions:
if(any(exclude))
  {
  dnonfitted <- data.frame(matrix(NA, nrow=length(curdnexclude), ncol=ncol(gof)))
  rownames(dnonfitted) <- curdnexclude
  colnames(dnonfitted) <- colnames(gof)
  gof <- rbind(gof, dnonfitted)
  }
}
# output:
dlf$gof <- gof
if(plot) distLgofPlot(dlf, quiet=quiet, ...)
dlf
} # end of function

