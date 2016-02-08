# Compute quantile of General Pareto Distribution fitted to sample
# Using the fExtremes package (Extreme Financial Market Data)
# Berry Boessenkool, Feb 2016, berry-b@gmx.de


#' GPD quantile of sample
#' 
#' Compute quantile of General Pareto Distribution fitted to sample by peak over treshold (POT) method
#' declaring treshold from truncation proportion using the fExtremes package (Extreme Financial Market Data)
#' 
#' @param x Vector with values. NAs are silently ignored
#' @param probs Probabilities of truncated (Peak over treshold) quantile
#' @param truncate Truncation percentage (proportion of sample discarded)
#' @param undertruncNA Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
#' @param method Method in \code{\link[fExtremes]{gpdFit}}, "mle" (maximum likelihood estimation) or "pwm" (probability-weighted moments). Note this is called type in the original function! DEFAULT: "pwm"
#' @param quiet Should messages be suppressed? DEFAULT: FALSE
#' @param \dots Further arguments passed to \code{\link[fExtremes]{gpdFit}}, except for \code{type}, as it is already specified by \code{method}.

#' @return Vector of quantile estimates for each value of probs
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @seealso \code{\link{distLquantile}}, \code{\link{q_evir}},
#' \code{\link[fExtremes]{gpdFit}}
#' @keywords distribution robust univar
#' @export
#' @importFrom fExtremes gpdFit
#' @importFrom berryFunctions pastec
#' @examples
#' 
#' library(fExtremes)
#' data(annMax)
#' p <- c(0.5, 0.8, 0.9, 0.99, 0.999)
#' d <- q_fExtremes(annMax, truncate=0,   p=p); d
#' d <- q_fExtremes(annMax, truncate=0.6, p=p); d
#' d <- q_fExtremes(annMax, truncate=0.6, p=p, under=FALSE); d
#' d <- q_fExtremes(annMax, truncate=1,   p=p); d
#' d <- q_fExtremes(annMax, truncate=1,   p=p, under=FALSE); d
#' 
q_fExtremes <- function(
x,
probs,
truncate,
undertruncNA=TRUE,
method="pwm",
quiet=FALSE,
...)
{
# Input control
x <- x[!is.na(x)]
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  if(!quiet) on.exit(message("Note in q_fExtremes: only first value of 'truncate' is used."), add=TRUE)
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
# prepare failure output:
failout <- rep(NA, length(probs))
names(failout) <- paste0(probs*100,"%")
# check probs for compliance with truncate:
if(all(probs < truncate) & !quiet)
   {
   on.exit(message("Note in q_fExtremes: 'probs' (", berryFunctions::pastec(probs),
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
# threshold from truncation proportion:
threshold <- if(truncate!=0) sort(x)[length(x)*truncate] else min(x)-1
# fitting: ---------------------------------------------------------------------
z <- try(z <- fExtremes::gpdFit(x, type=method, u=threshold, ...), silent=TRUE)
if(class(z)=="try-error")
  {
  if(!quiet) on.exit(message(
     "Note in q_fExtremes: fExtremes::gpdFit failed. Returning NAs. \n  Reason: ", z), add=TRUE)
  return(failout)
  }
output <- fExtremes::qgpd(p=probs2, xi=z@fit$par.ests["xi"], mu=z@parameter$u, beta=z@fit$par.ests["beta"])
output <- as.vector(output)
names(output) <- paste0(probs*100,"%")
# replace probs below truncation value with NA:
if(undertruncNA & any(probs < truncate) & !quiet)
  on.exit(message("Note in q_ismev: quantiles for probs below truncate (",
      berryFunctions::pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
output
}
