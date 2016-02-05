# Compute quantile of General Pareto Distribution fitted to sample
# Using the fExtremes package (Extreme Financial Market Data)
# Berry Boessenkool, Feb 2016, berry-b@gmx.de

q_fExtremes <- function(
x,           # Vector with values. NAs are silently ignored
probs,       # Probabilities of truncated (Peak over treshold) quantile
truncate,    # Truncation percentage (proportion of sample discarded)
undertruncNA=TRUE, # Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
method="pwm",# Method in \code{\link[fExtremes]{gpdFit}}, "mle" (maximum likelihood estimation) or "pwm" (probability-weighted moments). Note this is called type in the original function!
quiet=FALSE, # Should messages be suppressed?
...)         # Further arguments passed to \code{\link[fExtremes]{gpdFit}}, except for \code{type}, as it is already specified by \code{method}.
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
   on.exit(message("Note in q_fExtremes: 'probs' (",
        pastec(probs), ") must contain values that are larger than 'truncate' (",
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
z <- try(z <- fExtremes::gpdFit(x, type=method, u=threshold, ...))
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
                  pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
output
}
