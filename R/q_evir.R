# Compute quantile of General Pareto Distribution fitted to sample
# Using the established evir (extreme values in r) package
# Berry Boessenkool, End of July 2015
# berry-b@gmx.de


#' GPD quantile of sample
#' 
#' Compute quantile of General Pareto Distribution fitted to sample by peak over treshold (POT) method
#' declaring treshold from truncation proportion using the established evir (extreme values in r) package
#' 
#' @details \code{q_evir} opens and closes external device to suppress \code{evir::quant} plotting.\cr
#' \code{q_evir2} doesn't do that, but is a manual and slightly changed copy from the computing part of the code, thus
#' would not change in case \code{evir} is updated. It is also computationally much faster, as probs are vectorized.
#' 
#' @aliases q_evir q_evir2
#'
#' @param x Vector with values. NAs are silently ignored.
#' @param probs Probabilities of truncated (Peak over treshold) quantile.
#' @param truncate Truncation percentage (proportion of sample discarded).
#' @param undertruncNA Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
#' @param pngdev sink \code{evir::quant} graph output to file (is removed later) instead of openening \code{\link{dev.new}}, which also is closed later. Using TRUE avoids the graphics device popping up and disappearing again. Argument ignored in \code{q_evir2}. DEFAULT: TRUE
#' @param quantcat Show the cat messages of quant? Argument ignored in \code{q_evir2}. DEFAULT: FALSE
#' @param quiet Should messages be suppressed? DEFAULT: FALSE
#' @param quietgp Suppress q_evir gpd-optim failed notes? DEFAULT: quiet
#' @param method method in \code{\link{gpd}}, "ml" (maximum likelihood) or "pwm" (probability-weighted moments). Only used in \code{q_evir2}, ignored in \code{q_evir}. pwm yields higher results, but ml more often fails (gpd -> optim, Error in optim(theta, negloglik, hessian = TRUE, ..., tmp = excess) : non-finite finite-difference value [1]). DEFAULT: "pwm"
#' @param \dots Further arguments passed to \code{\link[evir]{quant}} or \code{\link[evir]{gpd}}

#' @return Vector of quantile estimates for each value of probs
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2015
#' @seealso \code{\link{distLquantile}}, \code{\link{q_ismev}}, \code{evir::\link[evir]{quant}}
#' @references \url{https://en.wikipedia.org/wiki/Pickands-Balkema-de_Haan_theorem} second
#'             extreme value theorem: tails of distributions tend towards GPD\cr
#'             \url{http://stats.stackexchange.com/questions/129438/different-quantiles-of-a-fitted-gpd-in-different-r-packages}
#' @keywords distribution robust univar
#' @export
#' @importFrom evir quant gpd
#' @importFrom berryFunctions pastec
#' @examples
#' 
#' library(lmomco)
#' library(berryFunctions)
#' library(evir)
#' 
#' # Data
#' set.seed(42); rnum <- rbeta(1000, 2, 7)*300
#' par(las=1)
#' evir::quant(rnum, p=0.95, models=70, end=1000)
#' axis(1, at=seq(-1000,0, length=6), labels=0:5/5, pos=par("usr")[4])
#' mtext("Proportion truncated", side=3, line=-3)
#' abline(v=-c(100,200), lty=2) 
#' # stable quantile estimate at the tail, but not for very high truncation
#' 
#' # To get one estimate, use q_evir or q_evir2 (faster, but code manually copypasted) 
#' q_evir(rnum, probs=0.95, truncate=0.8)
#' abline(h=q_evir(rnum, probs=0.95, truncate=0.8), lty=2)
#' 
#' evir::quant(rnum, p=0.999, models=70, end=1000)
#' # very high quantile estimates are much more uncertain, of course...
#' 
#' 
#' p <- c(0.5, 0.8, 0.9, 0.99, 0.999)
#' data(annMax)
#' q_evir(annMax, truncate=0,   pr=p)
#' q_evir(annMax, truncate=0.6, pr=p)
#' q_evir(annMax, truncate=0.6, pr=p, under=FALSE)
#' q_evir(annMax, truncate=1,   pr=p)
#' q_evir(annMax, truncate=1,   pr=p, under=FALSE)
#' 
#' library("ismev")
#' data(rain)
#' q_evir(rain, truncate=0,   pr=p)
#' q_evir(rain, truncate=0.6, pr=p)
#' q_evir(rain, truncate=0.6, pr=p, under=FALSE)
#' q_evir(rain, truncate=1,   pr=p)
#' q_evir(rain, truncate=1,   pr=p, under=FALSE)
#' 
#' q_evir(rain, truncate=0.2, pr=c(0.1,0.2,p))
#' q_evir2(rain, truncate=0.2, pr=c(0.1,0.2,p))
#' 
#' \dontrun{
#' ## Taken out from CRAN package check because it's slow
#' 
#' # I. GPD is not necessarily applicable for full samples!
#' d <- distLquantile(rnum, plot=TRUE, probs=0.99, nbest=18   )
#' # gpa (=GPD) performs badly here, compared to the other distributions, because:
#' distLgofPlot(distLfit(rnum, quiet=TRUE), ranks=FALSE) # GPA rmse is relatively high...
#' d <- distLquantile(rnum, plot=TRUE, probs=0.99, selection=c("wei", "gpa"),
#'               xlim=c(120, 270), ylim=c(0, 0.003), quiet=TRUE)
#' # and GPA has limited distribution support!! (dot at end)
#' # (computing quantiles or drawing random numbers will never exceed 184)
#' # Remember: the dist. parameters are in distLfit(rnum, sel="gpa")$param$gpa  or:
#' par <- lmomco::pargpa(lmomco::lmoms(rnum))
#' lmomco::quagpa(f=1, para=par)
#' 
#' 
#' # II. GPA performs better for distribution tails: 2nd extreme value theorem
#' # Theoretically, the tails of distributions converge to GPD
#' # see reference section for a link
#' distLquantile(rnum, plot=TRUE, probs=0.99, progbars=FALSE)
#' distLquantile(rnum, plot=TRUE, probs=0.99, truncate=0.8, progbars=FALSE)
#' distLquantile(rnum, plot=TRUE, probs=0.99, selection=c("wak", "gpa"),
#'               xlim=c(120, 270), ylim=c(0.6, 1), truncate=0.8, cdf=TRUE)
#' 
#' 
#' 
#' # III. evir::quant vs linear moments lmomco (implemented in extremeStat::distLquantile)
#' # a) nontruncated full sample:
#' comparequantiles <- function(trunc=0, xlim=c(0,210), ylim=c(trunc,1), ... )
#' {
#' probs <- seq(trunc,1, length=200)
#' d <- distLquantile(rnum, probs=probs, truncate=trunc, efast=TRUE, ...)
#' plot(ecdf(rnum), las=1, ylim=ylim, col=8, col.01line=NA, xlim=xlim, 
#'      xlab="Quantile estimate", main=paste(trunc*100, "% truncated"))
#' lines(d["quantileMean",], probs, col="orange")
#' lines(d["weighted1",], probs, col=1)
#' lines(d["gpa",],       probs, col=2)
#' lines(d["q_evir",],    probs, col=4)
#' legend("bottomright", c("ecdf","quantileMean","weighted","gpa","evir"), 
#'        col=c(8,"orange",1,2,4), lty=1)
#' invisible(d)
#' }
#' comparequantiles(ylim=lim0(1), xlim=c(0,210) )
#' abline(h=c(0.3,0.9), lty=2)
#' # evir GPD method underestimates 30% Quantile, overestimates at 90%
#' # lmomco gpa does the same, but not so extremely.
#' # Hardly surprising: GPD is not a beta distibution!
#' 
#' # Here are the numbers
#' dlq <- distLquantile(rnum, probs=c(0.3, 0.9), quiet=TRUE)
#' quantileMean(rnum, probs=c(0.3, 0.9))       # empirical
#' dlq["weighted1",]                           # fitted distributions
#' dlq["gpa",]                                 # GPD fitted by linear moments
#' q_evir(rnum, probs=c(0.3, 0.9), truncate=0) # GPD in package evir
#' 
#'          
#' # Now let's look at the tail of the beta-distributed random numbers:
#' # b) truncated Peak Over Treshold (POT) EV
#' comparequantiles(trunc=0.8, ylim=c(0.95, 1), xlim=c(140,210) )
#' abline(h=0.99, lty=2) 
#' # lmomco and evir yield the same result. Both slightly overestimate
#' 
#' comparequantiles(trunc=0.9, ylim=c(0.95, 1), xlim=c(140,210), quiet=TRUE )
#' d <- comparequantiles(trunc=0.98, xlim=c(155,255), quiet=TRUE )
#' 
#' 
#' # Comparison gpd methods:
#' set.seed(42); rnum <- rexp(500, rate=1/4)
#' samples <- lapply(1:100, function(n) sample(rnum,50)  )
#' results <- sapply(samples, q_evir2, probs=0.999, truncate=0.8)
#' sum(is.na(results))/100 # 19%
#' results2 <- sapply(samples, q_evir2, probs=0.999, truncate=0.8, method="pwm")
#' sum(is.na(results2))/100 # 0%
#' op <- par(mfrow=c(1,2))
#' plot(results, results2, las=1, xlab="ml", ylab="pwm", main="evir::gpd methods")
#' abline(a=0,b=1)
#' hist(results-results2, breaks=30, col="orange", main="hist differences")
#' # often, pwm is higher (negative diff), but sometimes ml can be MUCH higher 
#' 
#' # with different random numbers:
#' set.seed(42); rnum <- rbeta(1000, 2, 7)*300
#' samples <- lapply(1:100, function(n) sample(rnum,50)  )
#' results <- sapply(samples, q_evir2, probs=0.999, truncate=0.8)
#' sum(is.na(results))/100 # 45%
#' results2 <- sapply(samples, q_evir2, probs=0.999, truncate=0.8, method="pwm")
#' sum(is.na(results2))/100 # 0%
#' plot(results, results2, las=1, xlab="ml", ylab="pwm", main="evir::gpd methods")
#' abline(a=0,b=1)
#' hist(results-results2, breaks=30, col="orange", main="hist differences")
#' par(op)
#' }
#' 
q_evir <- function(
x,
probs,
truncate,
undertruncNA=TRUE,
pngdev=TRUE,
quantcat=FALSE,
quiet=FALSE,
quietgp=quiet,
method="pwm",
...)
{
# Input control
x <- x[!is.na(x)]
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  if(!quiet) on.exit(message("Note in q_evir: only first value of 'truncate' is used."), add=TRUE)
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
# failure output:
failout <- rep(NA, length(probs))
names(failout) <- paste0(probs*100,"%")
if(all(probs < truncate))
    {
    if(quiet) on.exit(message("Note in q_evir: 'probs' (",
    berryFunctions::pastec(probs), ") must contain values that are larger than 'truncate' (",
    truncate, "). Returning NAs."), add=TRUE)
    return(failout)
    }
# position of truncation (treshold)
pos <- length(x)*(1-truncate)
# quant always plots a graph, but we don't want it
tfile <- tempfile(fileext=".png")
if(pngdev) png(tfile) else dev.new()
at_end <- function(...) 
  {
  dummy <- dev.off();
  if(pngdev) unlink(tfile) else 
  if(dummy!=1) dev.set() # so active window is the one it used to be
  }
on.exit(at_end(), add=TRUE )
# actual computation with evir::quant
# object for quant cat results and message already given
allcats <- ""
OptimFailMessageGiven <- FALSE
if(pos<1)
  {
  if(!quiet) on.exit(message("Note in q_evir: Not enough values after truncation (",
                             pos,"). Returning NAs."), add=TRUE)
  return(failout)
  } else
  output <- sapply(probs, function(p)
  {
  cats <- capture.output(
    res <- try(evir::quant(x, p=p, start=pos, end=pos, models=1, ...), silent=TRUE)
    )
  assign("allcats", value=c(allcats, cats), envir=parent.env(environment()))
  if(class(res)=="try-error") 
  {
  if(!quietgp & !OptimFailMessageGiven)on.exit(message(
  "Note in q_evir: evir::quant-gpd-optim failed. Returning NAs.\n  Reason: ", res), add=TRUE)
  assign("OptimFailMessageGiven", value=TRUE, envir=parent.env(environment()))
  NA 
  } else
  as.numeric(res["qest",])
  })
# replace probs below truncation value with NA
if(undertruncNA & any(probs < truncate) & !quiet) 
  on.exit(message("Note in q_evir: quantiles for probs below truncate (",
      berryFunctions::pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Cat quantcats:
if(!quiet & quantcat) cat(allcats, sep="\n")
# Output result:
output
}



#' @export

# Copying only the computing part from evir::quant, Version: 1.7-3, Date: 2011-07-22
q_evir2 <- function(x, probs, truncate, undertruncNA=TRUE, 
                    pngdev=TRUE, quantcat=FALSE, quiet=FALSE, quietgp=quiet, method="pwm", ...)
  {
  x <- x[!is.na(x)]
  if(all(probs < truncate) & !quiet & undertruncNA) on.exit(message("Note in q_evir2: 'probs' (",
   berryFunctions::pastec(probs), ") must contain values that are larger than 'truncate' (",
    truncate, "). Returning NAs."), add=TRUE)
  
  pos <- length(x)*(1-truncate)
  if(pos<1)
  {
  if(!quiet) on.exit(message("Note in q_evir2: Not enough values after truncation (",
                             pos,"). Returning NAs."), add=TRUE)
  return(rep(NA, length(probs)))
  } else
  fit <- try(evir::gpd(x, nextremes=pos, method=method, ...), silent=TRUE)
  if(class(fit)=="try-error") 
  {
  if(!quietgp) on.exit(message("Note in q_evir2: evir::gpd-optim failed. Returning NAs.\n  Reason: ", fit), add=TRUE)
  return(rep(NA, length(probs)))
  } else
  lambda <- length(x)/fit$n.exceed
  a <- lambda * (1 - probs)
  gfunc <- function(a, xihat) (a^(-xihat) - 1)/xihat
  qest <- fit$threshold + fit$par.ests["beta"] * gfunc(a, fit$par.ests["xi"])
  # replace probs below truncation value with NA
  if(undertruncNA & any(probs < truncate) & !quiet) 
    on.exit(message("Note in q_evir2: quantiles for probs below truncate(",
          berryFunctions::pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
  if(undertruncNA) qest[probs < truncate] <- NA
  # Output
  unname(qest)
  }


# Wrong approach:
#q_evir3 <- function(x, probs, truncate, ...)
#  {
#  pos <- length(x)*(1-truncate)
#  fit <- evir::gpd(x, nextremes=pos, ...)
#  evir::qgpd(p=probs, xi=fit$par.ests["xi"], mu=0, beta=fit$par.ests["beta"])
#  }

