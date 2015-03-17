# Quantile of distributions fitted to a sample
# Berry Boessenkool, berry-b@gmx.de, March 2015

distLquantile <- function(
x=NULL,         # Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used.
probs=0:4/4,    # Numeric vector of probabilities with values in [0,1].
truncate=0,     # Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly.
type=NULL,      # Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used.
dlf=NULL,       # dlf object described in \code{\link{extremeStat}}. Ignored if x is given.
order=TRUE,     # Logical; if TRUE (the default) and length(type)>1, the output is ordered by dlf$gof, else by order of appearance in dlf$parameter
...             # Arguments passed to \code{\link{distLfit}} if x is given but dlf is not.
)
{
# input checks:
truncate <- truncate[1] # cannot be vectorized
if(truncate<0 | truncate>=1) stop("truncate must be a number between 0 and 1.")
if( is.null(x) & is.null(dlf)) stop("Either dlf or x must be given.")
if(!is.null(x) & truncate==0) dlf <- distLfit(dat=x, selection=type, ...)
probs2 <- probs
# truncation:
if(truncate!=0)
  {
  if(is.null(x)) x <- dlf$dat
  xtrunc <- sort(x)[ -1:-(truncate*length(x)) ]
  dlf <- distLfit(dat=xtrunc, selection=type, ...)
  probs2 <- (probs-truncate)/(1-truncate)
  probs2[probs < truncate] <- 0
  }
# available distributions
if(any(!type %in% names(dlf$parameter))) stop("Specified type is not available in dlf$parameter.")
available <- rownames(dlf$gof)
if(!is.null(type)) available <- available[available %in% type]
if(!is.null(type) & any(sort(available) != sort(type))) stop("Specified type is not available in dlf$gof.")
# Parameter object from list
param <- dlf$parameter
if(!is.null(type)) param <- param[type]
# distribution quantiles:
quan <- sapply(param, function(p) qlmomco(f=probs2, para=p))
if(length(probs2)==1) quan <- t(quan)
rownames(quan) <- paste0(probs*100,"%")
if(order & length(type)>1) quan <- quan[,available] # ordered by goodness of fit
quan[probs < truncate,] <- NA
quan
}
