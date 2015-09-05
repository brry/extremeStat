# Extreme value statistics for flood risk estimation
# Bootstrapping for confidence Interval
# Berry Boessenkool, 2015, berry-b@gmx.de

distLexBoot <- function(
dlf, # dlf object, as returned by \code{\link{distLextreme}}, is passed to \code{\link{distLextremePlot}}.
nbest=3, # Number of best fitted distribution functions in dlf for which bootstrapping is to be done. Overriden by \code{selection}.
selection=NULL, # Character vector with distribution function names to be used. Suggested to keep this low.
truncate=0, # Truncation of subsamples, see \code{\link{distLquantile}}
n=100, # Number of subsamples to be processed (computing time increases extraordinarily)
prop=0.8, # Proportion of sample to be used in each run
returnall=FALSE, # Return all simulations, instead of the aggregate confidence level?
conf.lev=0.95, # Confidence level (Proportion of subsamples within 'confidence interval'). Quantiles extracted from this value are passed to \code{\link[berryFunctions]{quantileMean}}
RPs=NULL, # Return Period vector, by default calculated internally based on log
plot=TRUE, # Plot results via \code{\link{distLextremePlot}}?
add=FALSE, # Add to existing plot?
log=TRUE, # Plot on a logarithmic axis.
progbars=TRUE, # Show progress bar for Monte Carlo simulation?
... # Further arguments passed to \code{\link{distLextremePlot}}
)
{
# Selection
if(is.null(selection)) selection <- rownames(dlf$gof)[1:nbest]
# Return period vector:
RPdef <- logSpaced(min=1, n=100, plot=FALSE, base=if(log) 1.1708 else 1, max=length(dlf$dat)*2)
if(is.null(RPs)) RPs <-  unique(round(RPdef, digits=2))
# subsample size:
sss <- round(length(dlf$dat)*prop)
# Actual computation for each subsample
if(progbars) replicate <- pbapply::pbreplicate
returnlev <- replicate(n=n,
   distLquantile(x=sample(dlf$dat, size=sss), selection=selection, truncate=truncate, 
   probs=1-1/RPs, empirical=FALSE, weighted=FALSE, trans=TRUE, progbars=FALSE, time=FALSE),
   simplify=FALSE)
# list restructuring
returnlev2 <- lapply(selection, function(i)
                               sapply(returnlev, "[", i=i, j=1:length(RPs))   )
names(returnlev2) <- selection
for(i in 1:length(returnlev2)) rownames(returnlev2[[i]]) <- RPs
# confidence Band calculation # toDo: allow vectorization
quant <- (1-conf.lev[1])/2
quant <- c(0+quant, 1-quant)
returnCI <- lapply(returnlev2, function(x)
                    apply(x, MARGIN=1, FUN=quantileMean, probs=quant, na.rm=TRUE))
# Plotting
if(plot & !add) distLextremePlot(dlf=dlf, selection=selection, nbest=nbest, log=log)
if(plot) for(i in length(returnCI):1)
ciBand(yu=returnCI[[i]][2,], yl=returnCI[[i]][1,], x=RPs, add=TRUE, colm=dlf$coldist[i], nastars=FALSE)
# Output
if(returnall) return(returnlev2) else return(returnCI)
}
