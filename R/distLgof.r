# distributions via linear moments
# Berry Boessenkool, Sept 2014
# calculate 'Goodness of Fit' measures in distLfit

#
distLgof <- function(
dlf, # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, datname, gofProp, parameter}
gofProp, # Overrides value in list. Proportion of highest values in \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
plot=TRUE, # Call \code{\link{distLgofplot}}?
progbars=TRUE, # Show progress bars for each loop?
ks=TRUE, # Include ks.test results in dlf$gof? Computing is much faster when FALSE
quiet=FALSE  # Should \code{\link{rmse}} warn about NA removal?
)
{
#browser()
# Progress bars
if( require(pbapply,quietly=TRUE) & progbars ) lapply <- pbapply::pblapply
# Objects from list:
dat <- dlf$dat
if(missing(gofProp)) gofProp <- dlf$gofProp
else if(length(gofProp)>1 | any(gofProp<0) | any(gofProp>1) ) stop("gofProp must be a single value between 0 and 1")
parameter <- dlf$parameter
dn <- names(parameter)
# Error check:
exclude <- sapply(parameter, function(x) if(!is.null(x)) any(is.na(x$para)) else TRUE)
if(any(exclude))
  {warning("The following distributions were excluded since no parameters were estimated:\n",
             paste(dn[exclude], collapse=", "))
  dn <- dn[!exclude]
  parameter <- parameter[!exclude] # not sure whether this is always good...
  }
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) print("performing ks.test:")
  ksA <- lapply(dn, function(d) ks.test(dat, paste0("cdf",d), parameter[[d]]) )
  ksP <- sapply(ksA, function(x) x$p.value   )
  ksD <- sapply(ksA, function(x) x$statistic )
  names(ksD) <- dn
  }
# CDFS for R2 on upper gofProp of data:
dat2 <- sort(dat, decreasing=TRUE)[  1:(gofProp*length(dat))  ]
if(progbars) print("calculating CDFs:")
tcdfs <- lapply(dn, function(d) plmomco(dat2,parameter[[d]]))
names(tcdfs) <- dn # Theoretical CumulatedDensityFunctions
ecdfs <- ecdf(dat)(dat2) # Empirical CDF
# Root Mean Square Error, R squared:
if( require(pbapply) & progbars ) sapply <- pbapply::pbsapply
if(progbars) print("calculating RMSE:")
RMSE <- sapply(dn, function(d)    rmse(tcdfs[[d]], ecdfs, quiet=quiet))
if(progbars) print("calculating R2:")
R2   <- sapply(dn, function(d) rsquare(tcdfs[[d]], ecdfs, quiet=TRUE))
# All into one data.frame:
gof <- data.frame(RMSE=RMSE, R2=R2)
if(ks) {gof$ksP=ksP; gof$ksD=ksD}
# order by GOF:
gof <- gof[ order(gof$RMSE), ]  # -gof$R2 # which measure should I sort by?
# output:
output <- list(dat=dat, datname=dlf$datname, gofProp=gofProp, parameter=parameter, gof=gof)
if(plot) distLgofplot(output)
output
} # end of function

