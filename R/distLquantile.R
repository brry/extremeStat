# Quantile of distributions fitted to a sample
# Berry Boessenkool, berry-b@gmx.de, March + July 2015

distLquantile <- function(
x=NULL,         # Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used.
probs=0:4/4,    # Numeric vector of probabilities with values in [0,1].
truncate=0,     # Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly.
selection=NULL, # Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used.
type=NULL,      # Overrides selection. Kept for backwards compatibility. (Changed to selection for internal consistency).
dlf=NULL,       # dlf object described in \code{\link{extremeStat}}. Ignored if x is given, or if truncate / selection do not match. Use this to save computing time for large datasets where you already have dlf.
order=TRUE,     # Sort results by GOF? If TRUE (the default) and length(selection)>1, the output is ordered by dlf$gof, else by order of appearance in selection (or dlf$parameter)
returndlf=FALSE,# Return full \code{dlf}list with output attached as element \code{quant}? If FALSE (the default), just the matrix with quantile estimates is returned.
plot=FALSE,     # Should \code{\link{distLplot}} be called?
cdf=FALSE,      # If TRUE, plot cumulated DF instead of probability density
lines=TRUE,     # Should vertical lines marking the quantiles be added?
linargs=NULL,   # Arguments passed to \code{\link{lines}}.
empirical=TRUE, # Add vertical line for empirical \code{\link{quantileMean}} (and include the result in the output matrix)?
ismev=empirical,# Include \code{evir::\link[ismev]{gpf.fit}} GPD quantile estimation? DEFAULT: empirical, so additional options can all be excluded with emp=F.
evir=empirical, # Include \code{evir::\link[evir]{quant}} GPD quantile estimation? Note that this temporarily creates a png image at the \code{getwd} if efast=FALSE. DEFAULT: empirical, so additional options can all be excluded with emp=F.
efast=TRUE,     # compute evir::quant only in the faster way (with  \code{\link{q_evir2})?
method=c("ml","pwm"), # Method in \code{\link{q_evir2}}, "ml" (maximum likelihood) or "pwm" (probability-weighted moments). Can also be both
weighted=empirical, # Include weighted averages across distribution functions to the output?
quiet=FALSE,    # Suppress notes?
quietss=quiet,  # Suppress sample size notes?
quietgp=quiet,  # Suppress q_evir gpd-optim failed notes?
...             # Arguments passed to \code{\link{distLfit}}.
)
{
# input checks: ----------------------------------------------------------------
internaldatname <- deparse(substitute(x))
truncate <- truncate[1] # cannot be vectorized
if(truncate<0 | truncate>=1) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
if( is.null(x) & is.null(dlf)) stop("Either dlf or x must be given.")
if(!is.null(type))
  {
  selection <- type
  if(!quiet) on.exit(message("Note in distLquantile: Argument 'type' overwrote 'selection'."), add=TRUE)
  }
if(!is.null(x)) if(is.list(x)) stop("x must be a vector. Possibly, you want to use dlf=",
                                    deparse(substitute(x)))
#
# Fit distribution functions to (truncated) sample: ----------------------------
if(!is.null(x) | any(dlf$truncate!=truncate) |  any(!selection %in% names(dlf$parameter)))
  {
  if(is.null(x)) 
    {
    x <- dlf$dat
    internaldatname <- dlf$datname
    if(!quiet)on.exit(message(
      "Note in distLquantile: selection or truncate given; distLfit is called from dlf$dat, as x is NULL."), add=TRUE)
    }
  # actual fitting:
  dlf <- distLfit(dat=x, datname=internaldatname, selection=selection, 
                  truncate=truncate, plot=plot, cdf=cdf, quiet=quiet, quietss=quietss, ...)
  }else
  {
  # reduce number of distfunctions analyzed if more were present in dlf argument:
  if(!is.null(selection))
    {
    dlf$parameter <- dlf$parameter[selection]
    dlf$gof <- dlf$gof[rownames(dlf$gof) %in% selection,]
    }
  }
#
# Empty output matrix: ---------------------------------------------------------
dn <- names(dlf$parameter)  # dn = distribution names
dn2 <- rownames(dlf$gof)
if(any(!dn2 %in% dn)) stop(pastec(dn2[!dn2 %in% dn]),
                             " available in dlf$gof, but not in dlf$parameter.")
if(any(!dn %in% dn2)) stop(pastec( dn[!dn %in% dn2]),
                             " available in dlf$parameter, but not in dlf$gof.")
if(order) dn <- dn2
output <- matrix(NA, ncol=length(probs), nrow=length(dn) )
colnames(output) <- paste0(probs*100,"%")
rownames(output) <- dn
#
# control for distributions that could not be fitted, e.g. kappa, but are in selection:
miss <- selection[!selection %in% dn]
if(length(miss)>0)
  {
  on.exit(message("Note in distLquantile: specified selection (", pastec(miss),
                  ") is not available in dlf$gof."), add=TRUE)
  m <- matrix(NA, nrow=length(miss), ncol=ncol(output))
  rownames(m) <- miss # always keep the same order if selection is given
  output <- cbind(output, m) # append missing distfuns at the end (which should be in correct place if order=F)
  }
#
# add rows for weighted averages of distribution functions
# and for GPD comparison methods
output <- rbind(output, quantileMean=NA, weighted1=NA, weighted2=NA, weighted3=NA,
            weightedc=NA, q_evir=NA, q_evir2_ml=NA, q_evir2_pwm=NA, q_ismev=NA)
#
# if input sample size is too small, return NA matrix:
if( length(dlf$dat)<5 )
  {
  if(!quietss) on.exit(message(
    "Note in distLquantile: sample size is too small to fit parameters (",
    length(dlf$dat),"). Returning NAs"), add=TRUE)
  if(returndlf) {dlf$quant <- output; return(dlf)}
  else return(output)
  }
#
# truncation probs update: -----------------------------------------------------
probs2 <- probs
if(truncate!=0)
  {
  if(all(probs < truncate) & !quiet) on.exit(message("Note in distLquantile: 'probs' (",
     pastec(probs),") must contain values that are larger than 'truncate (", 
          truncate, "). Returning NAs."), add=TRUE)
  probs2 <- (probs-truncate)/(1-truncate) # correct probabilities for smaller sample proportion
  probs2[probs < truncate] <- 0   # avoid negative values
  }
#
# distribution quantiles: ------------------------------------------------------
for(d in dn) if(!is.null(dlf$parameter[[d]])) {
               quantd <- lmomco::qlmomco(f=probs2, para=dlf$parameter[[d]])
               if(!is.null(quantd)) output[d,] <- quantd
               }
#
# Change results for probs below truncate to NA
if(truncate!=0) output[ , probs < truncate] <- NA
#
# add further quantile estimates -----------------------------------------------
# Empirical Quantiles:
if(empirical) output["quantileMean",] <- quantileMean(dlf$dat_full, probs=probs, truncate=truncate)
# evir::quant GPD estimates:
if(evir)  
  {
  if("ml" %in% method) output["q_evir2_ml",] <- q_evir2(  x=dlf$dat_full,
      probs=probs, truncate=truncate, method="ml", quiet=quiet, quietgp=quietgp)
  if("pwm" %in% method) output["q_evir2_pwm",] <- q_evir2(x=dlf$dat_full,
      probs=probs, truncate=truncate, method="pwm",quiet=quiet, quietgp=quietgp)
  if(!efast) output["q_evir",] <- q_evir(x=dlf$dat_full,
      probs=probs, truncate=truncate, quiet=quiet, quietgp=quietgp)
  }
# ismev::fit.gpf estimates:
if(ismev) output["q_ismev",] <- q_ismev(x=dlf$dat_full, probs=probs, truncate=truncate, quiet=quiet)
# Weighted quantile estimates:
if(weighted)
  {
  Qweighted <- function(Weightnr)
    {
    sapply(1:ncol(output), function(col_n)
      {
      vals <- output[rownames(dlf$gof), col_n]
      if(!any(is.finite(vals))) return(NA)
      weights <- dlf$gof[is.finite(vals),Weightnr]
      weights <- weights/sum(weights) # rescale to 1
      sum(vals[is.finite(vals)] * weights)
      })
    }
  output["weighted1",] <- Qweighted("weight1")
  output["weighted2",] <- Qweighted("weight2")
  output["weighted3",] <- Qweighted("weight3")
  output["weightedc",] <- Qweighted("weightc")
  }
#
# Plotting: Quantile lines: ----------------------------------------------------
if(plot & lines)
  {
  lfun <- if(cdf) plmomco else dlmomco
  use <- rownames(dlf$gof)[1:length(dlf$coldist)]
  for(i in length(use):1)
  {
  qval <- output[use[i],]
  qval <- qval[ is.finite(qval) ] # Inf ignored
  if(length(qval)>0) do.call(graphics::lines, args=owa(
                      list(x=qval, y=lfun(x=qval, para=dlf$parameter[[use[i]]]),
                  col=dlf$coldist[i], type="h"), linargs, "x","col","type"))
  }
  # empirical quantile added:
  if(empirical)
    {
    dd <- density(dlf$dat)
    qd <- output["quantileMean",]
    do.call(graphics::lines, args=owa(list(x=qd, y=approx(dd$x, dd$y, xout=qd)$y,
                                           type="h"), linargs, "x","y","type"))
    }
  }
# return output: ---------------------------------------------------------------
if(returndlf) {dlf$quant <- output; return(dlf)}
else return(output)
}

