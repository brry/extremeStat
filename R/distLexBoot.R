#' Bootstrapping uncertainty intervals for return periods
#' 
#' Calculates and plots bootstrap uncertainty intervals for \code{\link{plotLextreme}}.
#'
#' @details Has not been thoroughly tested yet. Bootstrapping defaults can probably be improved.
#' 
#' @return invisible dlf object, see \code{\link{printL}}.
#' Additional elements are: exBootSelection (names of distributions), 
#' exBootRPs (x values for plot)
#' exBootSim (all simulation results) and exBootCI (agregated into CI band). 
#' The last two are each a list with a matrix (RPs)
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2015 + Dec 2016
#' @seealso \code{\link{plotLexBoot}}, \code{\link{distLextreme}}
#' @keywords bootstrap montecarlo hplot dplot distribution ts
#' @export
#' @importFrom berryFunctions logSpaced quantileMean
#' 
#' @examples
#' 
#' data(annMax)
#' dlf <- distLextreme(annMax, selection=c("wak","gum","gev","nor"))
#' dlfB <- distLexBoot(dlf, nbest=4, conf.lev=0.5, n=10) # n low for quick example tests
#' plotLexBoot(dlfB)
#' 
#' @param dlf       \code{dlf} object, as returned by \code{\link{distLextreme}}, 
#'                  is passed to \code{\link{plotLextreme}}.
#' @param nbest     Number of best fitted distribution functions in dlf for which 
#'                  bootstrapping is to be done. Overriden by \code{selection}. DEFAULT: 3
#' @param selection Character vector with distribution function names to be used. 
#'                  Suggested to keep this low. DEFAULT: NULL
#' @param n         Number of subsamples to be processed 
#'                  (computing time increases extraordinarily). DEFAULT: 100
#' @param prop      Proportion of sample to be used in each run. DEFAULT: 0.8
#' @param conf.lev  Confidence level (Proportion of subsamples within 'confidence interval'). 
#'                  Quantiles extracted from this value are passed to 
#'                  \code{\link[berryFunctions]{quantileMean}}. DEFAULT: 0.95
#' @param RPs       Return Period vector, by default calculated internally based on 
#'                  value of \code{log}. DEFAULT: NULL
#' @param log       RPs suitable for plot on a logarithmic axis? DEFAULT: TRUE
#' @param progbars  Show progress bar for Monte Carlo simulation? DEFAULT: TRUE
#' @param \dots     Further arguments passed to \code{\link{distLquantile}}
#'                  like truncate, quiet=TRUE
#' 
distLexBoot <- function(
dlf,
nbest=3,
selection=NULL,
n=100,
prop=0.8,
conf.lev=0.95,
RPs=NULL,
log=TRUE,
progbars=TRUE,
...
)
{
# Selection
if(is.null(selection)) selection <- rownames(dlf$gof)[order(dlf$gof$RMSE)][1:nbest]
# Return period vector:
RPdef <- berryFunctions::logSpaced(min=1, n=100, plot=FALSE, 
                                   base=if(log) 1.1708 else 1, max=length(dlf$dat)*2)
if(is.null(RPs)) RPs <-  unique(round(RPdef, digits=2))
# subsample size:
sss <- round(length(dlf$dat)*prop)
# Actual computation for each subsample
if(progbars) replicate <- pbapply::pbreplicate
simQ_orig <- replicate(n=n,
   distLquantile(x=sample(dlf$dat, size=sss), selection=selection, order=FALSE,
   probs=1-1/RPs, empirical=FALSE, weighted=FALSE, progbars=FALSE, time=FALSE, ...),
   simplify=FALSE)
# list restructuring
simQ <- lapply(selection, function(i)
                               sapply(simQ_orig, "[", i=i, j=1:length(RPs))   )
names(simQ) <- selection
for(i in 1:length(simQ)) rownames(simQ[[i]]) <- RPs
# confidence Band calculation # toDo: allow vectorization
quant <- (1-conf.lev[1])/2
quant <- c(0+quant, 1-quant)
returnCI <- lapply(simQ, function(x)
       apply(x, MARGIN=1, FUN=berryFunctions::quantileMean, probs=quant, na.rm=TRUE))
# Output
dlf$exBootSim <- simQ
dlf$exBootCI <- returnCI
dlf$exBootRPs <- RPs
dlf$distnames <- selection
dlf$distcols <- berryFunctions::rainbow2(length(selection))
dlf$distselector <- "distLexBoot"
invisible(dlf)
}
