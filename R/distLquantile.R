# Quantile of distributions fitted to a sample
# Berry Boessenkool, berry-b@gmx.de, March 2015

distLquantile <- function(
x=NULL,         # Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used.
probs=0:4/4,    # Numeric vector of probabilities with values in [0,1].
truncate=0,     # Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly.
type=NULL,      # Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used.
dlf=NULL,       # dlf object described in \code{\link{extremeStat}}. Ignored if x is given.
order=TRUE,     # Sort results by GOF? If TRUE (the default) and length(type)>1, the output is ordered by dlf$gof, else by order of appearance in dlf$parameter
plot=FALSE,     # Should \code{\link{distLplot}} be called?
cdf=FALSE,      # If TRUE, plot cumulated DF instead of probability density
lines=TRUE,     # Should vertical lines marking the quantiles be added?
linargs=NULL,   # Arguments passed to \code{\link{lines}}.
empirical=TRUE, # Add vertical line for empirical \code{\link{quantileMean}}?
...             # Arguments passed to \code{\link{distLfit}} if x or truncate is given.
)
{
# input checks:
truncate <- truncate[1] # cannot be vectorized
if(truncate<0 | truncate>=1) stop("truncate must be a number between 0 and 1.")
if( is.null(x) & is.null(dlf)) stop("Either dlf or x must be given.")
if(!is.null(x) & truncate==0) dlf <- distLfit(dat=x, selection=type, plot=plot, cdf=cdf, ...)
probs2 <- probs
# truncation:
if(truncate!=0)
  {
  if(all(probs < truncate)) stop("probs must contain values that are bigger than truncate.")
  if(is.null(x)) x <- dlf$dat
  xtrunc <- sort(x)[ -1:-(truncate*length(x)) ]
  dlf <- distLfit(dat=xtrunc, selection=type, plot=plot, cdf=cdf, ...)
  probs2 <- (probs-truncate)/(1-truncate)
  probs2[probs < truncate] <- 0   # avoid negative values
  }
# available distributions
if(is.null(type)) type <- names(dlf$parameter)
tinp <- type %in% names(dlf$parameter)
if(any(!tinp)) message("note in distLquantile: Specified type (",
         paste(type[!tinp], collapse=", ") , ") is not available in dlf$parameter.")
type <- type[tinp]
available <- rownames(dlf$gof)
available <- available[available %in% type]
if(any(!type %in% available)) stop("Specified type is not available in dlf$gof.")
# Parameter object from list
param <- dlf$parameter[type]
if(length(type)==0) return() # stop("No type left")
# distribution quantiles:
quan <- sapply(param, function(p) qlmomco(f=probs2, para=p))
if(length(probs2)==1) quan <- t(quan)
rownames(quan) <- paste0(probs*100,"%")
# order by goodness of fit:
if(order & length(available)>1) quan <- quan[,available, drop=FALSE]
if(truncate!=0) quan[probs < truncate,] <- NA
# Quantile lines:
if(plot & lines)
  {
  lfun <- if(cdf) plmomco else dlmomco
  use <- available[1:length(dlf$coldist)]
  for(i in 1:length(use))
  do.call(graphics::lines, args=owa(list(x=quan[,use[i]], y=lfun(x=quan[,use[i]],
                                                  para=dlf$parameter[[use[i]]]),
                  col=dlf$coldist[i], type="h"), linargs, "x","y","col","type"))
  # empirical quantile added:
  if(empirical)
    {
    dd <- density(dlf$dat)
    qd <- quantileMean(dlf$dat, probs=probs, truncate=truncate)
    do.call(graphics::lines, args=owa(list(x=qd, y=approx(dd$x, dd$y, xout=qd)$y,
                                           type="h"), linargs, "x","y","type"))
    }
  }
# final output
quan
}
