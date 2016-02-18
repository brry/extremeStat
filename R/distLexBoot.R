#' Bootstrapping uncertainty intervals for return periods
#' 
#' Calculates and plots bootstrap uncertainty intervals for \code{\link{distLextremePlot}}.
#'
#' @details Has not been thoroughly tested yet. Bootstrapping defaults can probably be improved.
#' 
#' @return A list with (for each selection) a matrix with confidence intervals for RPs, or if returnall=TRUE, all the simulation results
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2015
#' @seealso \code{\link{distLextreme}}
#' @keywords bootstrap montecarlo hplot dplot distribution ts
#' @export
#' @importFrom berryFunctions logSpaced quantileMean ciBand
#' 
#' @examples
#' 
#' data(annMax)
#' dlf <- distLextreme(annMax, log=TRUE, selection=c("wak","gum","gev","nor"))
#' dleB <- distLexBoot(dlf, nbest=4, conf.lev=0.5)
#' 
#' @param dlf \code{dlf} object, as returned by \code{\link{distLextreme}}, is passed to \code{\link{distLextremePlot}}.
#' @param nbest Number of best fitted distribution functions in dlf for which bootstrapping is to be done. Overriden by \code{selection}. DEFAULT: 3
#' @param selection Character vector with distribution function names to be used. Suggested to keep this low. DEFAULT: NULL
#' @param truncate Truncation of subsamples, see \code{\link{distLquantile}}. DEFAULT: 0
#' @param n Number of subsamples to be processed (computing time increases extraordinarily). DEFAULT: 100
#' @param prop Proportion of sample to be used in each run. DEFAULT: 0.8
#' @param returnall Return all simulations, instead of the aggregate confidence level? DEFAULT: FALSE
#' @param conf.lev Confidence level (Proportion of subsamples within 'confidence interval'). Quantiles extracted from this value are passed to \code{\link[berryFunctions]{quantileMean}}. DEFAULT: 0.95
#' @param RPs Return Period vector, by default calculated internally based on log. DEFAULT: NULL
#' @param plot Plot results via \code{\link{distLextremePlot}}? DEFAULT: TRUE
#' @param add Add to existing plot? DEFAULT: FALSE
#' @param log Plot on a logarithmic axis. DEFAULT: TRUE
#' @param progbars Show progress bar for Monte Carlo simulation? DEFAULT: TRUE
#' @param \dots Further arguments passed to \code{\link{distLextremePlot}}
#' 
distLexBoot <- function(
dlf,
nbest=3,
selection=NULL,
truncate=0,
n=100,
prop=0.8,
returnall=FALSE,
conf.lev=0.95,
RPs=NULL,
plot=TRUE,
add=FALSE,
log=TRUE,
progbars=TRUE,
...
)
{
# Selection
if(is.null(selection)) selection <- rownames(dlf$gof)[1:nbest]
# Return period vector:
RPdef <- berryFunctions::logSpaced(min=1, n=100, plot=FALSE, base=if(log) 1.1708 else 1, max=length(dlf$dat)*2)
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
       apply(x, MARGIN=1, FUN=berryFunctions::quantileMean, probs=quant, na.rm=TRUE))
# Plotting
if(plot & !add) distLextremePlot(dlf=dlf, selection=selection, nbest=nbest, log=log)
if(plot) for(i in length(returnCI):1)
berryFunctions::ciBand(yu=returnCI[[i]][2,], yl=returnCI[[i]][1,], x=RPs, add=TRUE, colm=dlf$coldist[i], nastars=FALSE)
# Output
if(returnall) return(returnlev2) else return(returnCI)
}
