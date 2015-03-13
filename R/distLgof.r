# distributions via linear moments
# Berry Boessenkool, Sept 2014
# calculate 'Goodness of Fit' measures in distLfit

distLgof <- function(
dlf, # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, datname, gofProp, parameter}
gofProp, # Overrides value in list. Proportion of highest values in \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
plot=TRUE, # Call \code{\link{distLgofplot}}?
progbars=TRUE, # Show progress bars for each loop?
ks=TRUE, # Include ks.test results in dlf$gof? Computing is much faster when FALSE
quiet=FALSE  # Should \code{\link{rmse}} warn about NA removal?
)
{
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
  if(progbars) message("performing ks.test:")
  ksA <- lapply(dn, function(d) ks.test(dat, paste0("cdf",d), parameter[[d]]) )
  ksP <- sapply(ksA, function(x) x$p.value   )
  ksD <- sapply(ksA, function(x) x$statistic )
  names(ksD) <- dn
  }
# CDFS for R2 on upper gofProp of data:
dat2 <- sort(dat, decreasing=TRUE)[  1:(gofProp*length(dat))  ]
if(progbars) message("calculating CDFs:")
tcdfs <- lapply(dn, function(d) plmomco(dat2,parameter[[d]]))
names(tcdfs) <- dn # Theoretical CumulatedDensityFunctions
ecdfs <- ecdf(dat)(dat2) # Empirical CDF
# Root Mean Square Error, R squared:
if( require(pbapply) & progbars ) sapply <- pbapply::pbsapply
if(progbars) message("calculating RMSE:")
RMSE <- sapply(dn, function(d)    rmse(tcdfs[[d]], ecdfs, quiet=quiet))
if(progbars) message("calculating R2:")
R2   <- sapply(dn, function(d) rsquare(tcdfs[[d]], ecdfs, quiet=TRUE))
# All into one data.frame:
gof <- data.frame(RMSE=RMSE, R2=R2)
if(ks) {gof$ksP=ksP; gof$ksD=ksD}
# order by GOF:
gof <- gof[ order(gof$RMSE), ]  # -gof$R2 # which measure should I sort by?
# Weights for weighted average of return values:
# Weight including all distributions
gof$weight1 <- gof[,"RMSE"] # the lower, the better, the more weight
gof$weight1 <- max(gof$weight1)-gof$weight1+min(gof$weight1)  # with min or mean added,
gof$weight1 <- gof$weight1/sum(gof$weight1)  # the worst fit is not completely excluded
# Exclude worst fit
gof$weight2 <- gof[,"RMSE"]
gof$weight2 <- max(gof$weight2)-gof$weight2
gof$weight2 <- gof$weight2/sum(gof$weight2)
# use only best half:
gof$weight3 <- gof[,"RMSE"]
gof$weight3 <- max(gof$weight3)-gof$weight3
gof$weight3[(nrow(gof)/2):nrow(gof)] <- 0    # todo: test with selection length 1
gof$weight3 <- gof$weight3/sum(gof$weight3)
# output:
output <- list(dat=dat, datname=dlf$datname, gofProp=gofProp, parameter=parameter, gof=gof)
if(plot) distLgofplot(output)
output
} # end of function

