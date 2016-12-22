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
#' dlq <- distLquantile(annMax, weighted=FALSE, quiet=TRUE, probs=0.97, returnlist=TRUE)
#' dlq$quant
#' distLplot(dlq, qlines=TRUE) # per default best fitting distribution functions
#' distLplot(dlq, qlines=TRUE, qrow=c("wak","GPD*"), nbest=14)
#' #pdf("dummy.pdf", width=9)
#' distLplot(dlq, qlines=TRUE, qrow="GPD*", nbest=13, xlim=c(102,110), 
#'           qlinargs=list(lwd=3), qheights=seq(0.02, 0.005, len=14))
#' #dev.off()
#'
#'
#' \dontrun{
#' ## Taken out from CRAN package check because it's slow
#' 
#' # Sanity checks: important for very small samples:
#' x1 <- c(2.6, 2.5, 2.9, 3, 5, 2.7, 2.7, 5.7, 2.8, 3.1, 3.6, 2.6, 5.8, 5.6, 5.7, 5.3)
#' q1 <- distLquantile(x1, sanerange=c(0,500), sanevals=c(NA,500))
#' x2 <- c(6.1, 2.4, 4.1, 2.4, 6, 6.3, 2.9, 6.8, 3.5)
#' q2 <- distLquantile(x2, sanerange=c(0,500), sanevals=c(NA,500))
#' x3 <- c(4.4, 3, 1.8, 7.3, 2.1, 2.1, 1.8, 1.8)
#' q3 <- distLquantile(x3, sanerange=c(0,500), sanevals=c(NA,500))
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
#' trunc <- seq(0,0.1,len=200)
#' dd <- pbsapply(trunc, function(t) distLquantile(annMax, 
#'           selection="gpa", weight=FALSE, truncate=t, prob=0.99, quiet=T)[c(1,3),])
#   plot(trunc, dd[1,], type="o", las=1)
#' lines(trunc, dd[2,], type="o", col=2)
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
#' @param x Sample for which parametrical quantiles are to be calculated. 
#'          If it is NULL (the default), \code{dat} from \code{dlf} is used. DEFAULT: NULL
#' @param probs Numeric vector of probabilities with values in [0,1]. DEFAULT: c(0.8,0.9,0.99)
#' @param truncate Number between 0 and 1 (proportion of sample discarded). 
#'                 Censored quantile: fit to highest values only (truncate lower proportion of x). 
#'                 Probabilities are adjusted accordingly. DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, 
#'                  set this only via truncate, see Details of \code{\link{q_gpd}}. 
#'                  DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param weightc   Custom weights, see \code{\link{distLweights}}. DEFAULT: NA
#' @param sanerange Range outside of which results should be changed to \code{sanevals}.
#'                  This can capture numerical errors in small samples
#'                  (notably GPD_MLE_extRemes). If NA, this is ignored. Attention
#'                  if addinfo=T: the RMSE column is also checked and changed. DEFAULT: NA
#' @param sanevals  Values to be used below [1] and above [2] \code{sanerange}. DEFAULT: NA
#' @param selection Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. 
#'                  Can be a vector. If NULL (the default), all types present in 
#'                  dlf$parameter are used. DEFAULT: NULL
#' @param dlf dlf object described in \code{\link{extremeStat}}. Use this to save 
#'            computing time for large datasets where you already have dlf. DEFAULT: NULL
#' @param order Sort results by GOF? If FALSE, it is sorted by appearance in 
#'              selection (or dlf$parameter). DEFAULT: TRUE
#' @param returnlist Return full \code{dlf}list with output attached as element \code{quant}? 
#'                   If FALSE (the default), just the matrix with quantile estimates 
#'                   is returned. DEFAULT: FALSE
#' @param empirical Add empirical \code{\link{quantileMean}} in the output matrix 
#'                  and vertical lines? DEFAULT: TRUE
#' @param weighted Include weighted averages across distribution functions to the output?
#'                 DEFAULT: empirical, so additional options can all be excluded with emp=F.
#' @param gpd Include GPD quantile estimation via \code{\link{q_gpd}}? 
#'            Note that the 'GPD_LMO_lmomco' result differs slightly from 'gpa' if
#'            truncate=0. This comes from using x>threshold ('GPD_*') or
#'            x>=threshold ('gpa' and all other distributions in extremeStat). 
#'            DEFAULT: empirical
#' @param addinfo Should information about sample size and threshold be 
#'                \code{\link{rbind}ed} to the output? 
#'                Will also \code{\link{cbind}} RMSE values to the output.
#'                DEFAULT: FALSE
#' @param speed Compute \code{\link{q_gpd}} only for fast methods? 
#'              Don't accidentally set this to \code{FALSE} in simulations or 
#'              with large datasets! DEFAULT: TRUE
#' @param plot Should \code{\link{distLplot}} be called? DEFAULT: FALSE
#' @param plotargs List of arguments to be passed to \code{\link{distLplot}} 
#'                 like qlines, qheights, qrow, qlinargs, nbest, cdf, ...
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param ssquiet Suppress sample size notes? DEFAULT: quiet
#' @param ttquiet Suppress truncation!=threshold note? Note that \code{\link{q_gpd}} 
#'                is called with ttquiet=TRUE. DEFAULT: quiet
#' @param \dots Arguments passed to \code{\link{distLfit}}
#'              (and \code{\link{distLplot}} if plot=TRUE).
#'
distLquantile <- function(
x=NULL,
probs=c(0.8,0.9,0.99),
truncate=0,
threshold=berryFunctions::quantileMean(dlf$dat_full, truncate),
weightc=NA,
sanerange=NA,
sanevals=NA,
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
if(!is.null(dlf)) if(any(dlf$truncate!=truncate)|any(dlf$threshold!=threshold))
  {
  currentdlftruncate  <- dlf$truncate
  currentdlfthreshold <- dlf$threshold
  if(!quiet) on.exit(message("Note in distLquantile: truncate (",truncate,
       ") did not match dlf$truncate (",currentdlftruncate,
       "). Thresholds: ",toString(signif(c(threshold, currentdlfthreshold),7)),
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
                  truncate=truncate, threshold=threshold, order=order,
                  plot=FALSE, quiet=quiet, ssquiet=ssquiet, ...)
  }
# check selection
if(!is.null(selection))
  {
  if(any(!selection %in% names(dlf$parameter))) if(!quiet) on.exit(message(
   "Note in distLquantile: 'selection' (",selection[!selection %in% names(dlf$parameter)],
   ") is not in dlf$parameter. NAs will be returned for these distributions."), add=TRUE)

  # reduce number of distfunctions analyzed if more were present in dlf argument:
  dlf$parameter <- dlf$parameter[selection]
  dlf$gof <- dlf$gof[rownames(dlf$gof) %in% selection,]
  }
#
# check distribution names: ---------------------------------------------------------
dn <- rownames(dlf$gof)  # dn = distribution names
dn2 <- names(dlf$parameter)
if(any(!dn %in% dn2)) warning(toString(dn [!dn  %in% dn2]),
                             " available in dlf$gof, but not in dlf$parameter.")
if(any(!dn2%in% dn )) warning(toString(dn2[!dn2 %in% dn ]),
                             " available in dlf$parameter, but not in dlf$gof.")
# Empty output matrix: ---------------------------------------------------------
output <- matrix(NA, ncol=length(probs), nrow=length(dn) )
colnames(output) <- paste0(probs*100,"%")
rownames(output) <- dn
#
# control for distributions that could not be fitted, e.g. kappa, but are in selection:
miss <- selection[!selection %in% dn]
if(length(miss)>0)
  {
  on.exit(message("Note in distLquantile: specified selection (", toString(miss),
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
if(gpd) output <- rbind(output, 
    GPD_LMO_lmomco=NA, 
    GPD_LMO_extRemes=NA,
    GPD_PWM_evir=NA, 
    GPD_PWM_fExtremes=NA, 
    GPD_MLE_extRemes=NA,
    GPD_MLE_ismev=NA,
    GPD_MLE_evd=NA,
    GPD_MLE_Renext_Renouv=NA, 
    GPD_MLE_evir=NA, 
    GPD_MLE_fExtremes=NA,
    GPD_GML_extRemes=NA, 
    GPD_MLE_Renext_2par=NA,
    GPD_BAY_extRemes=NA)
if(addinfo) 
  {
  output <- cbind(output, RMSE=NA)
  output <- rbind(output, n_full=length(dlf$dat_full), n=length(dlf$dat), threshold=dlf$threshold)
  }
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
     toString(probs),") must contain values that are larger than 'truncate' (",
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
lenprob <- 1:length(probs)
for(d in dn) 
  {
  if(!is.finite(dlf$gof[d, "RMSE"])) next # maybe this is too restricitive...
  pard <- dlf$parameter[[d]]              # but it catches gno in skewed samples like
  if(is.null(pard)) next                  # c(2.4,2.7,2.3,2.5,2.2, 62.4 ,3.8,3.1) 
  if(inherits(pard, "try-error")) next
  if(d=="kap") if(pard$ifail!=0) next
  quantd <- try(lmomco::qlmomco(f=probs2, para=dlf$parameter[[d]]), silent=TRUE)
  if(inherits(quantd, "try-error")) next
  if(is.null(quantd)) next
  if(d=="kap") quantd[quantd<=pard$support[1] | quantd>=pard$support[2] ] <- NA
  output[d,lenprob] <- quantd
  }
#
# Change results for probs below truncate to NA
if(truncate!=0) output[ , probs < truncate] <- NA
# Add RMSE values if wanted:
if(addinfo) output[dn,"RMSE"] <- dlf$gof[dn,"RMSE"]
#
# add further quantile estimates -----------------------------------------------
# Empirical Quantiles:
if(empirical) output["quantileMean",lenprob] <- berryFunctions::quantileMean(dlf$dat_full,
                                             probs=probs, truncate=truncate)
# q_gpd estimates: -------------------------------------------------------------
if(gpd)
  {
  # inernal helper function:
  q_gpd_int <- function(pack, meth=NULL) q_gpd(package=pack, method=meth,
                        x=dlf$dat_full, probs=probs,
                        truncate=truncate, threshold=threshold, 
                        quiet=quiet, ttquiet=TRUE)
  #
  output["GPD_LMO_lmomco",lenprob]       <- q_gpd_int("lmomco")
  output["GPD_LMO_extRemes",lenprob]     <- q_gpd_int("extRemes", meth="Lmoments")
  output["GPD_PWM_evir",lenprob]         <- q_gpd_int("evir", meth="pwm")
  output["GPD_PWM_fExtremes",lenprob]    <- q_gpd_int("fExtremes", meth="pwm")
  output["GPD_MLE_extRemes",lenprob]     <- q_gpd_int("extRemes", meth="MLE")
  output["GPD_MLE_ismev",lenprob]        <- q_gpd_int("ismev")
  output["GPD_MLE_evd",lenprob]          <- q_gpd_int("evd")
  output["GPD_MLE_Renext_Renouv",lenprob]<- q_gpd_int("Renext", meth="r")
  output["GPD_MLE_evir",lenprob]         <- q_gpd_int("evir", meth="ml")
  output["GPD_MLE_fExtremes",lenprob]    <- q_gpd_int("fExtremes", meth="mle")
  output["GPD_GML_extRemes",lenprob]     <- q_gpd_int("extRemes", meth="GMLE")
  output["GPD_MLE_Renext_2par",lenprob]  <- q_gpd_int("Renext", meth="f")
if(!speed)output["GPD_BAY_extRemes",lenprob] <- q_gpd_int("extRemes", meth="Bayesian") # computes a while
  }
#
# sanity checks ----------------------------------------------------------------
if(!all(is.na(sanerange)))
  {
  sr <- range(sanerange, finite=TRUE)
  if(diff(sr)<sqrt(.Machine$double.eps)) stop("sanerange diff = 0. (",toString(sanerange),")")
  toosmall <- which(output<sr[1])
  toosmallnames <- rownames(which(output<sr[1], arr.ind=TRUE))
  toosmallnames <- toString(unique(toosmallnames))
  toosmallvals <- toString(output[toosmall])
  if(length(toosmall)>0)
  {
  if(!ssquiet) on.exit(message("Note in distLquantile: The following quantile ",
      "estimates are smaller than sanerange (",sr[1],") and changed to ",sanevals[1],":\n",
      toosmallnames, " with values: ", toosmallvals), add=TRUE)
  output[toosmall] <- sanevals[1]
  }
  toolarge <- which(output>sr[2])
  toolargenames <- rownames(which(output>sr[2], arr.ind=TRUE))
  toolargenames <- toString(unique(toolargenames))
  toolargevals <- toString(output[toolarge])
  if(length(toolarge)>0)
  {
  if(!ssquiet) on.exit(message("Note in distLquantile: The following quantile ",
      "estimates are larger than sanerange (",sr[2],") and changed to ",sanevals[2],":\n",
      toolargenames, " with values: ", toolargevals), add=TRUE)
  output[toolarge] <- sanevals[2]
  }
}
# Weighted quantile estimates: -------------------------------------------------
if(weighted) output <-  q_weighted(output, distLweights(dlf$gof, weightc=weightc))
#
dlf$quant <- output
# Plotting: Quantile lines: ----------------------------------------------------
if(plot) do.call(distLplot, berryFunctions::owa(c(list(dlf=dlf, qlines=TRUE),
                                                         list(...)), plotargs))
# return output: ---------------------------------------------------------------
if(returnlist) return(dlf) else return(output)
}
