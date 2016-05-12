#' distribution quantiles
#' 
#' Parametric quantiles of distributions fitted to a sample.
#' 
#' @details Very high quantiles (99\% and higher) need large sample sizes for
#' \code{\link{quantile}} to yield a robust estimate. Theoretically, at least
#' 1/(1-probs) values must be present, e.g. 10'000 for Q99.99\%. With smaller
#' sample sizes (eg n=35), they underestimate the actual (but unknown)
#' quantile. Parametric quantiles need only small sample sizes. They don't have
#' a systematical underestimation bias, but have higher variability.
#' 
#' @return Matrix with distribution quantile values (with NAs for probs below truncate), \cr
#' or, if returnlist=TRUE, a \code{dlf} list as described in \code{\link{extremeStat}}.
#' @note NAs are always removed from x in \code{\link{distLfit}}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, March + July 2015, Feb 2016
#' @seealso \code{\link{q_gpd}}, \code{\link{distLfit}}, Xian Zhou, Liuquan Sun and Haobo Ren (2000): Quantile estimation for left truncated and right censored data, Statistica Sinica 10
#'          \url{http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n411.pdf}\cr
#'          require("truncdist")
#' @references On GPD: \url{http://stats.stackexchange.com/questions/69438}
#' @keywords distribution robust univar
#' @export
#' @importFrom berryFunctions quantileMean owa
#' @importFrom lmomco plmomco dlmomco qlmomco
#' 
#' @examples
#'
#' data(annMax) # Annual Discharge Maxima (streamflow)
#'
#' distLquantile(annMax, emp=FALSE) # several distribution functions in lmomco
#' distLquantile(annMax, truncate=0.8, probs=0.95) # POT (annMax already block maxima)
#' distLquantile(annMax, probs=0.95, plot=TRUE, qlinargs=list(lwd=3), nbest=5, breaks=10)
#' # Parametric 95% quantile estimates range from 92 to 111!
#' # But the best fitting distributions all lie aroud 103.
#'
#' # compare General Pareto Fitting methods
#' # Theoretically, the tails of distributions converge to GPD (General Pareto)
#' # q_gpd compares several R packages for fitting and quantile estimation:
#' dlq <- distLquantile(annMax, weight=FALSE, quiet=TRUE, probs=0.97, returnlist=TRUE)
#' dlq$quant
#' distLplot(dlq, qlines=TRUE) # per default best fitting distribution functions
#' distLplot(dlq, qlines=TRUE, qrow=c("wak","q_gpd*"), nbest=14)
#' #pdf("dummy.pdf", width=9)
#' distLplot(dlq, qlines=TRUE, qrow="q_gpd*", nbest=13, xlim=c(102,110), 
#'           qlinargs=list(lwd=3), qheights=seq(0.02, 0.005, len=13))
#' #dev.off()
#'
#'
#' \dontrun{
#' ## Taken out from CRAN package check because it's slow
#' 
#' # weighted distribution quantiles are calculated by different weighting schemes:
#' dlf <- distLfit(annMax)
#' distLgofPlot(dlf, ranks=FALSE, weights=TRUE)
#' 
#' # If speed is important and parameters are already available, pass them via dlf:
#' distLquantile(dlf=dlf, probs=0:5/5, selection=c("wak","gev","kap"), order=FALSE)
#' distLquantile(dlf=dlf, truncate=0.3, returnlist=TRUE)$truncate
#'
#' # censored (truncated, trimmed) quantile, Peak Over Treshold (POT) method:
#' qwak <- distLquantile(annMax, sel="wak", prob=0.95, plot=TRUE, ylim=c(0,0.06), emp=FALSE)
#' qwak2 <-distLquantile(annMax, sel="wak", prob=0.95, truncate=0.6, plot=TRUE,
#'                      addinfo=FALSE, add=TRUE, coldist="blue", empirical=FALSE)
#'                      
#'
#' # Simulation of truncation effect
#' library(lmomco)
#' #set.seed(42)
#' rnum <- rlmomco(n=1e3, para=dlf$parameter$gev)
#' myprobs <- c(0.9, 0.95, 0.99, 0.999)
#' mytrunc <- seq(0, 0.9, length.out=20)
#' trunceffect <- sapply(mytrunc, function(mt) distLquantile(rnum, selection="gev",
#'                              probs=myprobs, truncate=mt, plot=FALSE, quiet=TRUE,
#'                              progbars=FALSE, empirical=FALSE)["gev",])
#' # If more values are truncated, the function runs faster
#'
#' op <- par(mfrow=c(2,1), mar=c(2,4.5,2,0.5), cex.main=1)
#' distLquantile(rnum, sel="gev", probs=myprobs, emp=FALSE, ylab="", xlab="", plot=TRUE)
#' distLquantile(rnum, sel="gev", probs=myprobs, emp=FALSE, addinfo=FALSE,
#'               truncate=0.3, add=TRUE, coldist=4, plot=TRUE)
#' legend("right", c("fitted GEV", "fitted with truncate=0.3"), lty=1, col=c(2,4),
#'        bg="white")
#' par(mar=c(3,4.5,3,0.5))
#' plot(mytrunc, trunceffect[1,], ylim=range(trunceffect), las=1, type="l",
#'      main=c("High quantiles of 1000 random numbers from gev distribution",
#'            "Estimation based on proportion of lower values truncated"),
#'      xlab="", ylab="parametrical quantile")
#' title(xlab="Proportion censored", mgp=c(1.8,1,0))
#' for(i in 2:4) lines(mytrunc, trunceffect[i,])
#' library("berryFunctions")
#' textField(rep(0.5,4), trunceffect[,11], paste0("Q",myprobs*100,"%") )
#' par(op)
#'
#'
#' set.seed(3); rnum <- rlmomco(n=1e3, para=dlf$parameter$gpa)
#' qd99 <- evir::quant(rnum, p=0.99, start=15, end=1000, ci=0.5, models=30)
#' axis(3, at=seq(-1000,0, length=6), labels=0:5/5, pos=par("usr")[3])
#' title(xlab="Proportion truncated", line=-3)
#' mytrunc <- seq(0, 0.9, length.out=30)
#' trunceffect <- sapply(mytrunc, function(mt) distLquantile(rnum, selection="gpa",
#'                       probs=0.99, truncate=mt, plot=FALSE, quiet=TRUE,
#'                       empirical=FALSE, gpd=TRUE))
#' lines(-1000*(1-mytrunc), trunceffect[1,], col=4)
#' lines(-1000*(1-mytrunc), trunceffect[2,], col=3) # interesting...
#' for(i in 3:13) lines(-1000*(1-mytrunc), trunceffect[i,], col=3) # interesting...
#'
#' # If you want the estimates only for one single truncation, use
#' q_gpd(rnum, probs=myprobs, truncate=0.5)
#' 
#' } # end dontrun
#' 
#' @param x Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used. DEFAULT: NULL
#' @param probs Numeric vector of probabilities with values in [0,1]. DEFAULT: c(0.8,0.9,0.99)
#' @param truncate Number between 0 and 1 (proportion of sample discarded). Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly. DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, set this only via truncate, see Details of \code{\link{q_gpd}}. DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param selection Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used. DEFAULT: NULL
#' @param dlf dlf object described in \code{\link{extremeStat}}. Use this to save computing time for large datasets where you already have dlf. DEFAULT: NULL
#' @param order Sort results by GOF? If TRUE (the default) and length(selection)>1, the output is ordered by dlf$gof, else by order of appearance in selection (or dlf$parameter). DEFAULT: TRUE
#' @param returnlist Return full \code{dlf}list with output attached as element \code{quant}? If FALSE (the default), just the matrix with quantile estimates is returned. DEFAULT: FALSE
#' @param empirical Add empirical \code{\link{quantileMean}} in the output matrix and vertical lines? DEFAULT: TRUE
#' @param weighted Include weighted averages across distribution functions to the output? DEFAULT: empirical, so additional options can all be excluded with emp=F.
#' @param gpd Include GPD quantile estimation via \code{\link{q_gpd}}? DEFAULT: empirical
#' @param addinfo Should information like sample size be \code{\link{rbind}ed} to the output? DEFAULT: FALSE
#' @param speed Compute \code{\link{q_gpd}} only for fast methods? Don't accidentally set this to \code{FALSE} in simulations or with large datasets! DEFAULT: TRUE
#' @param plot Should \code{\link{distLplot}} be called? DEFAULT: FALSE
#' @param plotargs List of arguments to be passed to \code{\link{distLplot}} like qlines, qheights, qrow, qlinargs, nbest, cdf, ...
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param ssquiet Suppress sample size notes? DEFAULT: quiet
#' @param ttquiet Suppress truncation!=threshold note? Note that \code{\link{q_gpd}} is called with ttquiet=TRUE. DEFAULT: quiet
#' @param \dots Arguments passed to \code{\link{distLfit}} (and \code{\link{distLplot}} if plot=TRUE).
#'
distLquantile <- function(
x=NULL,
probs=c(0.8,0.9,0.99),
truncate=0,
threshold=berryFunctions::quantileMean(dlf$dat_full, truncate),
selection=NULL,
dlf=NULL,
order=TRUE,
returnlist=FALSE,
empirical=TRUE,
weighted=empirical,
gpd=empirical,
addinfo=FALSE,
speed=TRUE,
plot=FALSE,
plotargs=NULL,
quiet=FALSE,
ssquiet=quiet,
ttquiet=quiet,
...
)
{
# input checks: ----------------------------------------------------------------
internaldatname <- deparse(substitute(x))
truncate <- truncate[1] # cannot be vectorized
if(truncate<0 | truncate>=1) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
if( is.null(x) &  is.null(dlf)) stop("Either dlf or x must be given.")
if(!is.null(x) & !is.null(dlf)) stop("Either dlf or x must be given, but not both.")
if(!is.null(x)) if(is.list(x)) stop("x must be a vector. Possibly, you want to use dlf=",
                                    deparse(substitute(x)))
#
# Fit distribution functions to (truncated) sample: ----------------------------
# check truncate
if(!is.null(dlf)) if(any(dlf$truncate!=truncate)|any(any(dlf$threshold!=threshold)))
  {
  currentdlftruncate  <- dlf$truncate
  currentdlfthreshold <- dlf$threshold
  if(!quiet) on.exit(message("Note in distLquantile: truncate (",truncate,
       ") did not match dlf$truncate (",currentdlftruncate,
       "). Thresholds: ",pastec(signif(c(threshold, currentdlfthreshold),7)),
       ".\n   distLfit is called with the original dlf$dat."), add=TRUE)
  x <- dlf$dat
  internaldatname <- dlf$datname
  dlf <- NULL
  }
# actual fitting:
if(is.null(dlf))
  {
  # threshold initialization impossible if dlf = NULL
  if(is.na(threshold)) threshold <- berryFunctions::quantileMean(x, truncate)
  dlf <- distLfit(dat=x, datname=internaldatname, selection=selection,
                  truncate=truncate, threshold=threshold,
                  plot=FALSE, quiet=quiet, ssquiet=ssquiet, ...)
  }
# check selection
if(any(!selection %in% names(dlf$parameter))) if(!quiet) on.exit(message(
   "Note in distLquantile: 'selection' (",selection[!selection %in% names(dlf$parameter)],
   ") is not in dlf$parameter. NAs will be returned for these distributions."), add=TRUE)

# reduce number of distfunctions analyzed if more were present in dlf argument:
if(!is.null(selection))
  {
  dlf$parameter <- dlf$parameter[selection]
  dlf$gof <- dlf$gof[rownames(dlf$gof) %in% selection,]
  }
#
# Empty output matrix: ---------------------------------------------------------
dn <- names(dlf$parameter)  # dn = distribution names
dn2 <- rownames(dlf$gof)
if(any(!dn2 %in% dn)) warning(pastec(dn2[!dn2 %in% dn]),
                             " available in dlf$gof, but not in dlf$parameter.")
if(any(!dn %in% dn2)) warning(pastec( dn[!dn %in% dn2]),
                             " available in dlf$parameter, but not in dlf$gof.")
if(order) dn <- dn2
output <- matrix(NA, ncol=length(probs), nrow=length(dn) )
colnames(output) <- paste0(probs*100,"%")
rownames(output) <- dn
#
# control for distributions that could not be fitted, e.g. kappa, but are in selection:
miss <- selection[!selection %in% dn]
if(length(miss)>0)
  {
  on.exit(message("Note in distLquantile: specified selection (", pastec(miss),
                  ") is not available in dlf$gof."), add=TRUE)
  m <- matrix(NA, nrow=length(miss), ncol=ncol(output))
  rownames(m) <- miss # always keep the same order if selection is given
  output <- rbind(output, m) # append missing distfuns at the end (which should be in correct place if order=F)
  }
#
# add rows for weighted averages of distribution functions
# and for GPD comparison methods
if(empirical) output <- rbind(output, quantileMean=NA)
if(weighted) output <- rbind(output,weighted1=NA, weighted2=NA, weighted3=NA, weightedc=NA)
if(gpd) output <- rbind(output, q_gpd_evir_pwm=NA, q_gpd_evir_ml=NA,
    q_gpd_evd=NA, q_gpd_extRemes_MLE=NA, q_gpd_extRemes_GMLE=NA,
    q_gpd_extRemes_Bayesian=NA, q_gpd_extRemes_Lmoments=NA,
    q_gpd_fExtremes_pwm=NA, q_gpd_fExtremes_mle=NA, q_gpd_ismev=NA,
    q_gpd_Renext_r=NA, q_gpd_Renext_f=NA)
if(addinfo) output <- rbind(output, n_full=length(dlf$dat_full), n=length(dlf$dat), threshold=dlf$threshold)
#
# if input sample size is too small, return NA matrix:
if( length(dlf$dat)<5 )
  {
  if(!ssquiet) on.exit(message(
    "Note in distLquantile: sample size is too small to fit parameters (",
    length(dlf$dat),"). Returning NAs"), add=TRUE)
  if(returnlist) {dlf$quant <- output; return(dlf)}
  else return(output)
  }
#
# truncation probs update: -----------------------------------------------------
probs2 <- probs
if(truncate!=0)
  {
  if(all(probs < truncate) & !quiet) on.exit(message("Note in distLquantile: 'probs' (",
     pastec(probs),") must contain values that are larger than 'truncate' (", 
          truncate, "). Returning NAs."), add=TRUE)
  probs2 <- (probs-truncate)/(1-truncate) # correct probabilities for smaller sample proportion
  probs2[probs < truncate] <- 0   # avoid negative values
  }
#
# Threshold check:
normalthr <- berryFunctions::quantileMean(dlf$dat_full, truncate)
if(threshold != normalthr)
    {
    probs2 <- probs
    if(!ttquiet) on.exit(message("Note in distLquantile: threshold (",threshold,
    ") is not equal to threshold computed from truncate (",normalthr,
    ").\n  Probabilities are not corrected for truncation!"), add=TRUE)
    }
#
# distribution quantiles: ------------------------------------------------------
# This is the actual work...
for(d in dn) if(!is.null(dlf$parameter[[d]])) {
               quantd <- lmomco::qlmomco(f=probs2, para=dlf$parameter[[d]])
               if(!is.null(quantd)) output[d,] <- quantd
               }
#
# Change results for probs below truncate to NA
if(truncate!=0) output[ , probs < truncate] <- NA
#
# add further quantile estimates -----------------------------------------------
# Empirical Quantiles:
if(empirical) output["quantileMean",] <- berryFunctions::quantileMean(dlf$dat_full,
                                             probs=probs, truncate=truncate)
# Weighted quantile estimates:
if(weighted)
  {
  Qweighted <- function(Weightnr)
    {
    sapply(1:ncol(output), function(col_n)
      {
      vals <- output[rownames(dlf$gof), col_n]
      if(!any(is.finite(vals))) return(NA)
      weights <- dlf$gof[is.finite(vals),Weightnr]
      weights <- weights/sum(weights) # rescale to 1
      sum(vals[is.finite(vals)] * weights)
      })
    }
  output["weighted1",] <- Qweighted("weight1")
  output["weighted2",] <- Qweighted("weight2")
  output["weighted3",] <- Qweighted("weight3")
  output["weightedc",] <- Qweighted("weightc")
  }
# q_gpd estimates:
if(gpd)
  {
  # inernal helper function:
  q_gpd_int <- function(pack, meth=NULL) q_gpd(package=pack, method=meth,
                        x=dlf$dat_full, probs=probs,
                        truncate=truncate, threshold=threshold, 
                        quiet=quiet, ttquiet=TRUE)
  #
  output["q_gpd_evir_pwm",]      <- q_gpd_int("evir", meth="pwm")
  output["q_gpd_evir_ml",]       <- q_gpd_int("evir", meth="ml")
  output["q_gpd_evd",]           <- q_gpd_int("evd")
  output["q_gpd_extRemes_MLE",]  <- q_gpd_int("extRemes", meth="MLE")
  output["q_gpd_extRemes_GMLE",] <- q_gpd_int("extRemes", meth="GMLE")
if(!speed) output["q_gpd_extRemes_Bayesian",] <- q_gpd_int("extRemes", meth="Bayesian") # computes a while
  output["q_gpd_extRemes_Lmoments",] <- q_gpd_int("extRemes", meth="Lmoments")
  output["q_gpd_fExtremes_pwm",] <- q_gpd_int("fExtremes", meth="pwm")
  output["q_gpd_fExtremes_mle",] <- q_gpd_int("fExtremes", meth="mle")
  output["q_gpd_ismev",]         <- q_gpd_int("ismev")
  output["q_gpd_Renext_r",]      <- q_gpd_int("Renext", meth="r")
  output["q_gpd_Renext_f",]      <- q_gpd_int("Renext", meth="f")
  }
#
dlf$quant <- output
# Plotting: Quantile lines: ----------------------------------------------------
if(plot) do.call(distLplot, berryFunctions::owa(c(list(dlf=dlf, qlines=TRUE), list(...)), plotargs))
# return output: ---------------------------------------------------------------
if(returnlist) return(dlf) else return(output)
}
