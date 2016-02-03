# Compute quantile of General Pareto Distribution fitted to sample
# Using the ismev package (Introduction to Statistical Modeling of Extreme Values)
# Berry Boessenkool, Feb 2016, berry-b@gmx.de
# gpdq is not exported by ismev, so this is an ugly hack

q_ismev <- function(
x,           # Vector with values. NAs are silently ignored
probs,       # Probabilities of truncated (Peak over treshold) quantile
truncate,    # Truncation percentage (proportion of sample discarded)
returnlist=FALSE, # Return result from \code{\link[ismev]{gpd.fit}} with the quantiles added to the list as element \code{quant}
undertruncNA=TRUE, # Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
quiet=FALSE, # Should messages be suppressed?
...)         # Further arguments passed to \code{\link[ismev]{gpd.fit}}
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
   on.exit(message("Note in q_ismev: 'probs' (",
        pastec(probs), ") must contain values that are larger than 'truncate' (",
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
z <- try(gpd.fit(x, threshold=threshold, show=FALSE))
if(class(z)=="try-error")
  {
  if(!quiet) on.exit(message(
     "Note in q_ismev: ismev::gpd failed. Returning NAs. \n  Reason: ", z), add=TRUE)
  if(returnlist) return(list(z="not fitted", quant=failout)) else return(failout)
  }
output <- ismev:::gpdq(a=z$mle, u=z$threshold, p=1-probs2)
names(output) <- paste0(probs*100,"%")
# replace probs below truncation value with NA:
if(undertruncNA & any(probs < truncate) & !quiet)
  on.exit(message("Note in q_ismev: quantiles for probs below truncate (",
                  pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
if(returnlist) output <- list(z, quant=output)
output
}
