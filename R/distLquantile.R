# Quantile of distributions fitted to a sample
# Berry Boessenkool, berry-b@gmx.de, March + July 2015

distLquantile <- function(
x=NULL,         # Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used.
probs=0:4/4,    # Numeric vector of probabilities with values in [0,1].
truncate=0,     # Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly.
selection=NULL, # Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used.
type=NULL,      # Overrides selection. Kept for backwards compatibility. (Changed to selection for internal consistency).
dlf=NULL,       # dlf object described in \code{\link{extremeStat}}. Ignored if x is given, or if truncate / selection do not match. Use this to save computing time for large datasets where you already have dlf.
order=TRUE,     # Sort results by GOF? If TRUE (the default) and length(type)>1, the output is ordered by dlf$gof, else by order of appearance in dlf$parameter
plot=FALSE,     # Should \code{\link{distLplot}} be called?
cdf=FALSE,      # If TRUE, plot cumulated DF instead of probability density
lines=TRUE,     # Should vertical lines marking the quantiles be added?
linargs=NULL,   # Arguments passed to \code{\link{lines}}.
empirical=TRUE, # Add vertical line for empirical \code{\link{quantileMean}} (and include the result in the output matrix)?
evir=empirical, # Include \code{evir::\link[evir]{quant}} GPD quantile estimation? Note that this temporarily creates a png image at the \code{getwd}. DEFAULT: empirical, so additional options can all be excluded with emp=F.
efast=FALSE,    # compute evir::quant in a faster way (with  \code{\link{q_evir2})?
weighted=empirical, # Include weighted averages across distribution functions to the output?
quiet=FALSE,    # Suppress notes?
trans=quiet,    # Suppress note about transposing? # Option and message will be removed around the end of 2015.
...             # Arguments passed to \code{\link{distLfit}}.
)
{
# input checks: ----------------------------------------------------------------
internaldatname <- deparse(substitute(x))
# temporary warning:
if(!trans) message("Please note: distLquantile output has been transposed since Version 0.4.23 from 2015-07-18!")
truncate <- truncate[1] # cannot be vectorized
if(truncate<0 | truncate>=1) stop("truncate must be a number between 0 and 1.")
if( is.null(x) & is.null(dlf)) stop("Either dlf or x must be given.")
if(!is.null(type))
  {
  selection <- type
  if(!quiet) message("note in distLquantile: Argument 'type' overwrote 'selection'.")
  }
if(!is.null(x)) if(is.list(x)) stop("x must be a vector. Possibly, you want to use dlf =", deparse(substitute(x)))
#
# Fit distribution functions to (truncated) sample: ----------------------------
if(!is.null(x) | any(dlf$truncate!=truncate) |  any(!selection %in% names(dlf$parameter)))
  {
  dlf <- distLfit(dat=x, datname=internaldatname, selection=selection, 
                  truncate=truncate, plot=plot, cdf=cdf, quiet=quiet, ...)
  }else
  {
  # reduce number of distfunctions analyzed if more were present in dlf argument:
  dlf$parameter <- dlf$parameter[selection]
  dlf$gof <- dlf$gof[rownames(dlf$gof) %in% selection,]
  }
#
# Empty output matrix: ---------------------------------------------------------
output <- matrix(NA, nrow=length(probs), ncol=length(dlf$parameter) )
rownames(output) <- paste0(probs*100,"%")
colnames(output) <- names(dlf$parameter)
#
# control for distributions that could not be fitted, e.g. kappa:
miss <- selection[!selection %in% rownames(dlf$gof)]
miss <- miss[!miss %in% colnames(output)]
if(length(miss)>0)
  {
  message("note in distLquantile: specified selection (", pastec(miss),") is not available in dlf$gof.")
  m <- matrix(NA, ncol=length(miss), nrow=nrow(output))
  output <- cbind(output, m)
  colnames(output) <- selection # always keep the same order if selection is given
  }
#
# if input sample size is too small, return NA matrix:
if( length(dlf$dat)<5 )
  {
  if(!quiet) message("note in distLquantile: sample size is too small to fit parameters. Returning NAs")
  #if(empirical) output <- cbind(output, quantileMean=quantileMean(dlf$dat, probs=probs))
  if(empirical) output <- cbind(output, quantileMean=NA)
  return(t(output))
  }
#
# truncation probs update: -----------------------------------------------------
probs2 <- probs
if(truncate!=0)
  {
  if(all(probs < truncate) & !quiet) message("Note in distLquantile: 'probs' (",
     pastec(probs),") must contain values that are larger than 'truncate (", 
          truncate, "). Returning NAs.")
  probs2 <- (probs-truncate)/(1-truncate) # correct probabilities for smaller sample proportion
  probs2[probs < truncate] <- 0   # avoid negative values
  }
#
# distribution quantiles: ------------------------------------------------------
dn <- colnames(output)
for(d in dn) if(!is.null(dlf$parameter[[d]])) output[,d] <-
                                   qlmomco(f=probs2, para=dlf$parameter[[d]])
#
# Optional operations: ---------------------------------------------------------
# Change results for probs below truncate to NA
if(truncate!=0) output[probs < truncate,] <- NA
#
# order by goodness of fit:
if(order) output <- output[,rownames(dlf$gof), drop=FALSE]
# append missing distfuns at the end (which should be in correct place if order=F)
if(order & length(miss)>0)
  {
  m <- matrix(NA, ncol=length(miss), nrow=nrow(output))
  colnames(m) <- miss
  output <- cbind(output, m)
  }
#
# add further quantile estimates -----------------------------------------------
# Empirical Quantiles:
if(empirical) output <- cbind(output, quantileMean=quantileMean(
                                  dlf$dat_full, probs=probs, truncate=truncate))
# evir::quant GPD estimates:
if(evir)  
  {
  q_evir_sel <- if(efast)  q_evir2  else  q_evir
  output <-cbind(output, q_evir=q_evir_sel(
                   x=dlf$dat_full, probs=probs, truncate=truncate, quiet=quiet))
  }
# Weighted quantile estimates
if(missing(weighted) & !is.null(selection)) # false if not specified explicitely but
  if(length(selection) < 4) weighted <- FALSE     # only few distfuns are chosen 
if(weighted)
  {
  Qweighted <- function(Weightnr)
    {
    sapply(1:nrow(output), function(row_n)
      {
      vals <- output[row_n,rownames(dlf$gof)]
      if(!any(is.finite(vals))) return(NA)
      weights <- dlf$gof[is.finite(vals),Weightnr]
      weights <- weights/sum(weights) # rescale to 1
      sum(vals[is.finite(vals)] * weights)
      })
    }
  output <- cbind(output, weighted1=Qweighted("weight1"))
  output <- cbind(output, weighted2=Qweighted("weight2"))
  output <- cbind(output, weighted3=Qweighted("weight3"))
  }
#
# Plotting: Quantile lines: ----------------------------------------------------
if(plot & lines)
  {
  lfun <- if(cdf) plmomco else dlmomco
  use <- rownames(dlf$gof)[1:length(dlf$coldist)]
  for(i in length(use):1)
  {
  qval <- output[,use[i]]
  qval <- qval[ is.finite(qval) ] # Inf ignored
  if(length(qval)>0) do.call(graphics::lines, args=owa(
                      list(x=qval, y=lfun(x=qval, para=dlf$parameter[[use[i]]]),
                  col=dlf$coldist[i], type="h"), linargs, "x","y","col","type"))
  }
  # empirical quantile added:
  if(empirical)
    {
    dd <- density(dlf$dat)
    qd <- output[,"quantileMean"]
    do.call(graphics::lines, args=owa(list(x=qd, y=approx(dd$x, dd$y, xout=qd)$y,
                                           type="h"), linargs, "x","y","type"))
    }
  }
# return output: ---------------------------------------------------------------
t(output)
}

