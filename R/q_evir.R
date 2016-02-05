# Compute quantile of General Pareto Distribution fitted to sample
# Using the established evir (extreme values in r) package
# Berry Boessenkool, End of July 2015
# berry-b@gmx.de

q_evir <- function(
x,           # Vector with values. NAs are silently ignored.
probs,       # Probabilities of truncated (Peak over treshold) quantile.
truncate,    # Truncation percentage (proportion of sample discarded).
undertruncNA=TRUE, # Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
pngdev=TRUE, # sink \code{evir::quant} graph output to file (is removed later) instead of openening \code{\link{dev.new}}, which also is closed later. Using TRUE avoids the graphics device showing briefly.
quantcat=FALSE, # Show the cat messages of quant?
quiet=FALSE, # Should messages be suppressed?
quietgp=quiet,  # Suppress q_evir gpd-optim failed notes?
method="pwm", # method in \code{\link{gpd}}, "ml" (maximum likelihood) or "pwm" (probability-weighted moments). Only used in \code{q_evir2}, ignored in \code{q_evir}. pwm yields higher results, but ml more often fails (gpd -> optim, Error in optim(theta, negloglik, hessian = TRUE, ..., tmp = excess) : non-finite finite-difference value [1])
...)         # Further arguments passed to \code{\link[evir]{quant}} or \code{\link[evir]{gpd}}
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
    pastec(probs), ") must contain values that are larger than 'truncate' (", 
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
    res <- try(quant(x, p=p, start=pos, end=pos, models=1, ...), silent=TRUE)
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
                  pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Cat quantcats:
if(!quiet & quantcat) cat(allcats, sep="\n")
# Output result:
output
}





# Copying only the computing part from evir::quant, Version: 1.7-3, Date: 2011-07-22
q_evir2 <- function(x, probs, truncate, undertruncNA=TRUE, 
                    pngdev=TRUE, quantcat=FALSE, quiet=FALSE, quietgp=quiet, method="pwm", ...)
  {
  x <- x[!is.na(x)]
  if(all(probs < truncate) & !quiet & undertruncNA) on.exit(message("Note in q_evir2: 'probs' (",
    pastec(probs), ") must contain values that are larger than 'truncate' (", 
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
                    pastec(probs[probs<truncate]),") replaced with NAs."), add=TRUE)
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

