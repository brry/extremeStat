# Compute quantile of General Pareto Distribution fitted to sample
# Using the established evir (extreme values in r) package
# Berry Boessenkool, End of July 2015
# berry-b@gmx.de

q_evir <- function(
x,           # Vector with values
probs,       # Probabilities of truncated (Peak over treshold) quantile
truncate,    # Truncation percentage (proportion of sample discarded)
pngdev=TRUE, # sink \code{evir::quant} graph output to file (is removed later) instead of openening \code{\link{dev.new}}, which also is closed later. Using TRUE avoids the graphics device showing, but \code{getwd()} must be writable.
...)         # Further arguments passed to \code{\link[evir]{quant}}
{
# Input control
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  warning("only one value used for 'truncate'.")
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
# position of truncation (treshold)
pos <- length(x)*(1-truncate)
# quant always plots a graph, but we don't want it
if(pngdev) png("q_evir_dummy_plot.png") else dev.new()
# actual computation with evir::quant
output <- sapply(probs, function(p)
  {
  res <- quant(x, p=p, start=pos, end=pos, models=1, ...)
  as.numeric(res["qest",])
  })
dummy <- dev.off()
if(pngdev) unlink("q_evir_dummy_plot.png") else if(dummy!=1) dev.set() # so active window is the one it used to be
# Output result:
output
}


# Copying only the computing part from evir::quant, Version: 1.7-3, Date: 2011-07-22
q_evir2 <- function(x, probs, truncate, ...)
  {
  pos <- length(x)*(1-truncate)
  fit <- evir::gpd(x, nextremes=pos, ...)
  lambda <- length(x)/fit$n.exceed
  a <- lambda * (1 - probs)
  gfunc <- function(a, xihat) (a^(-xihat) - 1)/xihat
  qest <- fit$threshold + fit$par.ests["beta"] * gfunc(a, fit$par.ests["xi"])
  qest
  }


# Wrong approach:
#q_evir3 <- function(x, probs, truncate, ...)
#  {
#  pos <- length(x)*(1-truncate)
#  fit <- evir::gpd(x, nextremes=pos, ...)
#  evir::qgpd(p=probs, xi=fit$par.ests["xi"], mu=0, beta=fit$par.ests["beta"])
#  }

