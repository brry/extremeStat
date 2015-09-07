# distributions via linear moments
# Berry Boessenkool, Sept 2014, July 2015
# calculate 'Goodness of Fit' measures in distLfit

distLgof <- function(
dlf, # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, datname, gofProp, parameter}
gofProp, # Overrides value in list. Proportion of highest values in \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
plot=TRUE, # Call \code{\link{distLgofPlot}}?
progbars=length(dlf$dat)>200, # Show progress bars for each loop?
ks=TRUE, # Include ks.test results in dlf$gof? Computing is much faster when FALSE
quiet=FALSE, # Suppress notes?
... # Further arguments passed to \code{\link{distLgofPlot}}
)
{
# Progress bars
if(quiet) progbars <- FALSE
if(progbars) lapply <- pbapply::pblapply
# Objects from list:
dat <- dlf$dat
if(missing(gofProp)) gofProp <- dlf$gofProp
else if(length(gofProp)>1 | any(gofProp<0) | any(gofProp>1) ) stop("gofProp must be a single value between 0 and 1.")
parameter <- dlf$parameter
dn <- names(parameter)
# Error check:
exclude <- sapply(parameter, function(x) 
  {
  if(is.null(x)) return(TRUE)
  # if("ifail" %in% names(x)) if(x$ifail != 0) return(TRUE) ## restriction too tight
  if(is.null(plmomco(mean(dat),x))) return(TRUE)  # *)
  any(is.na(x$para))
  })                    #  *): CDF cannot be computed for kappa in Dresden example
if(any(exclude))
  {
  curdnexclude <- dn[exclude]
  if(!quiet) on.exit(message("Note in distLgof: The following distributions were excluded since no parameters were estimated:\n",
             pastec(curdnexclude)))
  dn <- dn[!exclude]
  # parameter <- parameter[!exclude] # not sure whether this is always good...
}
if(length(dn)<1) on.exit(message("No fitted distributions in dlf."))
if(length(dn)<2&!quiet) on.exit(message("Note in distLgof: Only ", pastec(dn), 
                                " was fitted, thus GOF can't be compared."))
if(ks)
  {
  # Kolmogorov-Smirnov test:
  if(progbars) message("Performing ks.test:")
  ksA <- lapply(dn, function(d) ks.test(dat, paste0("cdf",d), parameter[[d]]) )
  ksP <- sapply(ksA, function(x) x$p.value   )
  ksD <- sapply(ksA, function(x) x$statistic )
  names(ksD) <- dn
  }
# CDFS for R2 on upper gofProp of data:
dat2 <- sort(dat, decreasing=TRUE)[  1:(gofProp*length(dat))  ]
if(progbars) message("Calculating CDFs:")
tcdfs <- lapply(dn, function(d) plmomco(dat2,parameter[[d]]))
names(tcdfs) <- dn # Theoretical CumulatedDensityFunctions
ecdfs <- ecdf(dat)(dat2) # Empirical CDF
# Root Mean Square Error, R squared:
if(progbars) sapply <- pbapply::pbsapply
if(progbars) message("Calculating RMSE:")
RMSE <- sapply(dn, function(d)    rmse(tcdfs[[d]], ecdfs, quiet=TRUE))
if(progbars) message("Calculating R2:")
R2   <- sapply(dn, function(d) rsquare(tcdfs[[d]], ecdfs, quiet=TRUE))
if(!quiet)
  {
  nNA <- base::sapply(tcdfs, function(x) sum(is.na(x)))
  if(any(nNA>0)) 
    {
    dNA <- paste(paste0(dn[nNA>0], " (", nNA[nNA>0], ")"), collapse=", ")
    on.exit(message("Note in distLgof: NAs removed in CDF (limited support region?): ", 
                    dNA, " of ", length(tcdfs[[1]]), " values."))
    }
  }
# All into one data.frame:
gof <- data.frame(RMSE=RMSE, R2=R2)
if(!all(dim(gof) == 0)) # dim = 0,0 if all distributions are excluded
{
if(ks) {gof$ksP=ksP; gof$ksD=ksD}
# order by GOF:
gof <- gof[ order(gof$RMSE), ]  # -gof$R2 # which measure should I sort by?
# Weights for weighted average of return values:
# Weight including all distributions
gof$weight1 <- gof[,"RMSE"] # the lower, the better, the more weight
gof$weight1 <- max(gof$weight1)-gof$weight1+min(gof$weight1)  # with min or mean added,
gof$weight1 <- gof$weight1/sum(gof$weight1)  # the worst fit is not completely excluded
# Exclude worst fit (needs 2 or more distributions in selection to work)
gof$weight2 <- gof[,"RMSE"]
gof$weight2 <- max(gof$weight2)-gof$weight2
gof$weight2 <- gof$weight2/sum(gof$weight2)
# use only best half (needs 4 or more dists - technically, 3 should be enough...)
gof$weight3 <- gof[,"RMSE"]
gof$weight3 <- max(gof$weight3)-gof$weight3
gof$weight3[(nrow(gof)/2):nrow(gof)] <- 0
gof$weight3 <- gof$weight3/sum(gof$weight3)
} else
{
gof <- data.frame(matrix(NA, ncol=5, nrow=length(dlf$parameter) ))
colnames(gof) <- c("RMSE","R2","weight1","weight2","weight3")
rownames(gof) <-  names(dlf$parameter)
}
# output:
output <- list(dat=dat, datname=dlf$datname, gofProp=gofProp, parameter=parameter, gof=gof)
if(plot) distLgofPlot(output, quiet=quiet, ...)
output
} # end of function

