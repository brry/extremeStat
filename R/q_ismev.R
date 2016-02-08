# Compute quantile of General Pareto Distribution fitted to sample
# Using the ismev package (Introduction to Statistical Modeling of Extreme Values)
# Berry Boessenkool, Feb 2016, berry-b@gmx.de
# gpdq is not exported by ismev, so this is an ugly hack


#' GPD quantile of sample
#' 
#' Compute quantile of General Pareto Distribution fitted to sample by peak over treshold (POT) method
#' declaring treshold from truncation proportion using the ismev package (Introduction to Statistical Modeling of Extreme Values)
#' 
#' @param x Vector with values. NAs are silently ignored
#' @param probs Probabilities of truncated (Peak over treshold) quantile
#' @param truncate Truncation percentage (proportion of sample discarded)
#' @param returnlist Return result from \code{\link[ismev]{gpd.fit}} with the quantiles added to the list as element \code{quant}. DEFAULT: FALSE
#' @param undertruncNA Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
#' @param quiet Should messages be suppressed? DEFAULT: FALSE
#' @param \dots Further arguments passed to \code{\link[ismev]{gpd.fit}}

#' @return if(returnlist): gpd.fit list with element \code{quant} added \cr
#' else: Vector of quantile estimates for each value of probs
#' @note gpdq is not exported by ismev, so this is an ugly hack
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @seealso \code{\link{distLquantile}}, \code{\link{q_evir}},
#' \code{\link[ismev]{gpd.fit}}
#' @keywords distribution robust univar
#' @export
#' @importFrom ismev gpd.fit
#' @importFrom berryFunctions pastec
#' @examples
#' 
#' library(ismev)
#' data(rain)
#' data(annMax)
#' p <- c(0.5, 0.8, 0.9, 0.99, 0.999)
#' d <- q_ismev(annMax, truncate=0,   p=p); d
#' d <- q_ismev(annMax, truncate=0,   p=p, returnlist=TRUE); d
#' d <- q_ismev(annMax, truncate=0.6, p=p); d
#' d <- q_ismev(annMax, truncate=0.6, p=p, under=FALSE); d
#' d <- q_ismev(annMax, truncate=0.6, p=p, returnlist=TRUE); d
#' d <- q_ismev(annMax, truncate=1,   p=p); d
#' d <- q_ismev(annMax, truncate=1,   p=p, returnlist=TRUE); d
#' d <- q_ismev(annMax, truncate=1,   p=p, under=FALSE); d
#' 
#' 
#' d <- q_ismev(rain, truncate=0,   p=p); d
#' d <- q_ismev(rain, truncate=0,   p=p, returnlist=TRUE); tail(d,1)
#' d <- q_ismev(rain, truncate=0.6, p=p); d
#' d <- q_ismev(rain, truncate=0.6, p=p, under=FALSE); d
#' d <- q_ismev(rain, truncate=0.6, p=p, returnlist=TRUE); tail(d,1)
#' d <- q_ismev(rain, truncate=1,   p=p); d
#' d <- q_ismev(rain, truncate=1,   p=p, returnlist=TRUE); tail(d,1)
#' d <- q_ismev(rain, truncate=1,   p=p, under=FALSE); d
#' 
q_ismev <- function(
x,
probs,
truncate,
returnlist=FALSE,
undertruncNA=TRUE,
quiet=FALSE,
...)
{
# Input control
x <- x[!is.na(x)]
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  if(!quiet) on.exit(message("Note in q_ismev: only first value of 'truncate' is used."), add=TRUE)
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
# prepare failure output:
failout <- rep(NA, length(probs))
names(failout) <- paste0(probs*100,"%")
# check probs for compliance with truncate:
if(all(probs < truncate) & !quiet)
   {
   on.exit(message("Note in q_ismev: 'probs' (", berryFunctions::pastec(probs),
        ") must contain values that are larger than 'truncate' (",
        truncate, "). Returning NAs."), add=TRUE)
   if(returnlist) return(list(z="not fitted", quant=failout)) else return(failout)
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
z <- try(ismev::gpd.fit(x, threshold=threshold, show=FALSE), silent=TRUE)
if(class(z)=="try-error")
  {
  if(!quiet) on.exit(message(
     "Note in q_ismev: ismev::gpd.fit failed. Returning NAs. \n  Reason: ", z), add=TRUE)
  if(returnlist) return(list(z="not fitted", quant=failout)) else return(failout)
  }
output <- ismev:::gpdq(a=z$mle, u=z$threshold, p=1-probs2)
names(output) <- paste0(probs*100,"%")
# replace probs below truncation value with NA:
if(undertruncNA & any(probs < truncate) & !quiet)
  on.exit(message("Note in q_ismev: quantiles for probs below truncate (",
     berryFunctions::pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
if(returnlist) output <- list(z, quant=output)
output
}
