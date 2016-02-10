# Compute quantile of General Pareto Distribution fitted to sample
# Using one of several possible packages
# Berry Boessenkool, Feb 2016, berry-b@gmx.de


#' GPD quantile of sample
#' 
#' Compute quantile of General Pareto Distribution fitted to sample by peak over treshold (POT) method
#' using treshold from truncation proportion,
#' comparing several R packages doing this
#'
#' @details Depending on the value of "package", this fits the GPD using \cr
#' \code{evir::\link[evir]{gpd}}\cr
#' \code{extRemes::\link[extRemes]{fevd}}\cr
#' \code{fExtremes::\link[fExtremes]{gpdFit}}\cr
#' \code{ismev::\link[ismev]{gpd.fit}}\cr
#' The \code{method} defaults (and other possibilities) are \cr
#' evir: "pwm" (probability-weighted moments), or "ml" (maximum likelihood) \cr
#' extRemes: "MLE", or "GMLE", "Bayesian", "Lmoments" \cr
#' fExtremes: "pwm", or "mle"\cr
#' ismev: none, only Maximum-likelihood fitting implemented \cr
#' The Quantiles are always given with \code{probs} in regard to the full (uncensored) sample.
#' If e.g. truncate is 0.90, the distribution function is fitted to the top 10\% of the sample.
#' The 95th percentile of the full sample is equivalent to the 50\% quantile of the subsample acutally used for fitting.
#' For computation, the probabilities are internally updated with \code{p2=(p-t)/(1-t)}
#' but labelled with the original \code{p}.
#' If you truncate 90\% of the sample, you cannot compute the 70th percentile anymore,
#' thus \code{undertruncNA} should be left to TRUE. \cr
#' If not exported by the packages, the quantile functions are extracted from their current (Feb 2016) source code.
#' 
#' @param x Vector with numeric values. NAs are silently ignored.
#' @param probs Probabilities of truncated (Peak over treshold) quantile. DEFAULT: c(0.8,0.9,0.99)
#' @param truncate Truncation percentage (proportion of sample discarded). DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, set this only via truncate, see Details. DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param package Character string naming package to be used. One of c("evir","extRemes","fExtremes","ismev"). DEFAULT: "extRemes"
#' @param method \code{method} passed to the fitting function, if applicable. Defaults are internally specified, depending on \code{package}, if left to the DEFAULT: NULL.
#' @param returnlist Return result from the fitting funtion (if provided) with the quantiles added to the list as element \code{quant}. DEFAULT: FALSE
#' @param undertruncNA Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
#' @param quiet Should messages from this function be suppressed? DEFAULT: FALSE
#' @param suppresswarnings Should warnings be suppressed via \code{\link{options}(warn=-1)}? The usual type of warning is: NAs produced in log(...). DEFAULT: quiet
#' @param \dots Further arguments passed to the fitting funtion listed in section Details.

#' @return if(returnlist): list with element \code{quant} added \cr
#' else: Named vector of quantile estimates for each value of \code{probs}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @seealso \code{\link{distLquantile}} which compares results for all packages
#' @keywords distribution robust univar
#' @export
#' @importFrom evir quant gpd
#' @importFrom extRemes fevd qevd
#' @importFrom fExtremes gpdFit
#' @importFrom ismev gpd.fit
#' @importFrom berryFunctions pastec quantileMean
#'
#' @examples
#'  # q_gpd examples are currently in external file!
#'  stop("Berry, put the q_gpd examples back here!")
#'
q_gpd <- function(
x,
probs=c(0.8,0.9,0.99),
truncate=0,
threshold=berryFunctions::quantileMean(x, truncate),
package="extRemes",
method=NULL,
returnlist=FALSE,
undertruncNA=TRUE,
quiet=FALSE,
suppresswarnings=quiet,
...)
{
# Input control: ---------------------------------------------------------------
if(length(package)!=1) stop("package must have length 1, not ", length(package))
if(!package %in% c("evir","extRemes","fExtremes","ismev")) stop("package (",
     package, ") must be one of c('evir','extRemes','fExtremes','ismev')")
x <- x[!is.na(x)]
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  if(!quiet) on.exit(message("Note in q_gpd: only first value of 'truncate' is used."), add=TRUE)
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
# prepare failure output:
failout <- rep(NA, length(probs))
names(failout) <- paste0(probs*100,"%")
if(returnlist) failout <- c(list(z="not fitted"), list(quant=failout))
# check probs for compliance with truncate:
if(all(probs < truncate) & !quiet)
   {
   on.exit(message("Note in q_gpd: 'probs' (", berryFunctions::pastec(probs),
        ") must contain values that are larger than 'truncate' (",
        truncate, "). Returning NAs."), add=TRUE)
   return(failout)
   }
# truncation probs update: -----------------------------------------------------
probs2 <- probs
if(truncate!=0)
  {
  probs2 <- (probs-truncate)/(1-truncate) # correct probabilities for smaller sample proportion
  probs2[probs < truncate] <- 0   # avoid negative values
  }
#
# threshold warning:
if(!missing(threshold))
  {
  normalthr <- berryFunctions::quantileMean(x, truncate)
  if(threshold != normalthr)
    {
    probs2 <- probs
    if(!quiet) on.exit(message("Note in q_gpd: threshold (",threshold,
    ") is not equal to threshold computed from truncate (",normalthr,
    ").\n  Probabilities are not corrected for truncation!"), add=TRUE)
    }
  }
#
#
# fitting: ---------------------------------------------------------------------
if(suppresswarnings) {oo <- options(warn=-1); on.exit(options(oo), add=TRUE)}
# function to perform in case of failure. Yields useful error message and returns NAs:
failfun <- function(z, fitfun) {
  if(!quiet) on.exit(message(
     "Note in q_gpd: ",fitfun," failed. Returning NAs. \n  Reason: ", z), add=TRUE)
  return(failout)
  }
# actual fitting:
if(package=="evir") ##################
{
if(missing(method)) method <- "pwm"
pos <- sum(x > threshold)
z <- try(evir::gpd(x, nextremes=pos, method=method, ...), silent=TRUE)
if(class(z)=="try-error") return(failfun(z, "evir::gpd"))
} else
if(package=="extRemes") ##################
{
if(missing(method)) method <- "MLE"
z <- try(extRemes::fevd(x, method=method, type="GP", threshold=threshold, ...), silent=TRUE)
if(class(z)=="try-error") return(failfun(z, "extRemes::fevd"))
} else
if(package=="fExtremes") ##################
{
if(missing(method)) method <- "pwm"
z <- try(z <- fExtremes::gpdFit(x, type=method, u=threshold, ...), silent=TRUE)
if(class(z)=="try-error") return(failfun(z, "fExtremes::gpdFit"))
} else
if(package=="ismev") ##################
{
z <- try(ismev::gpd.fit(x, threshold=threshold, show=FALSE, ...), silent=TRUE)
if(class(z)=="try-error") return(failfun(z, "ismev::gpd.fit"))
} else
stop("package ", package, "is not in the options. This is a bug. Please report to berry-b@gmx.de.")
#
#
# quantile computing: ----------------------------------------------------------
if(package=="evir") ##################
{
  # computing part from evir::quant, Version: 1.7-3, Date: 2011-07-22
  lambda <- length(x)/z$n.exceed
  a <- lambda * (1 - probs)
  gfunc <- function(a, xihat) (a^(-xihat) - 1)/xihat
  output <- z$threshold + z$par.ests["beta"] * gfunc(a, z$par.ests["xi"])
} else
if(package=="extRemes") ##################
{
  # Get parameters from result:
  if(z$method=="Bayesian")
  {
  scale <- exp(mean(z$results[,"log.scale"]))
  shape <-     mean(z$results[,"shape"])
  } else
  if(z$method=="Lmoments")
  {
  scale <- z$results["scale"]
  shape <- z$results["shape"]
  } else
  {
  scale <- z$results$par["scale"]
  shape <- z$results$par["shape"]
  }
  probs2[probs2==0] <- NA
  output <- extRemes::qevd(p=probs2, scale=scale, shape=shape, threshold=z$threshold, type="GP")
} else
if(package=="fExtremes") ##################
{
  output <- fExtremes::qgpd(p=probs2, xi=z@fit$par.ests["xi"], mu=z@parameter$u, beta=z@fit$par.ests["beta"])
  output <- as.vector(output)
  z2 <- lapply(slotNames(z), getElement, object=z)
  names(z2) <- slotNames(z)
  z2$BerryWarning <- "transformed into list from Formal class 'fGPDFIT' [package 'fExtremes'] with 8 slots"
  z <- z2
} else
if(package=="ismev") ##################
{
  # from ismev Version 1.40 Date 2009-14-07, Published 2014-12-24    ismev:::gpdq
  ismev_gpdq <- function(a,u,p) u + (a[1] * (p^(-a[2]) - 1))/a[2]
  output <- ismev_gpdq(a=z$mle, u=z$threshold, p=1-probs2)
} else
stop("package ", package, "is not in the options. This is a bug. Please report to berry-b@gmx.de.")
#
#
# output formatting, checks: ---------------------------------------------------
names(output) <- paste0(probs*100,"%")
# replace probs below truncation value with NA:
if(undertruncNA & any(probs < truncate) & !quiet)
  on.exit(message("Note in q_gpd: quantiles for probs (",
     berryFunctions::pastec(probs[probs<=truncate]),
     ") below truncate (",truncate,") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
if(returnlist) output <- c(z, list(quant=output))
output
}

