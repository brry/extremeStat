#' Bootstrapping uncertainty intervals for return periods
#' 
#' Calculates and plots bootstrap uncertainty intervals for \code{\link{plotLextreme}}.
#' 
#' @details Has not been thoroughly tested yet. Bootstrapping defaults can probably be improved.
#' 
#' @return invisible dlf object, see \code{\link{printL}}.
#' Additional elements are: exBootCL (confidence level),
#' exBootRPs (x values for plot)
#' exBootSim (all simulation results) and exBootCI (aggregated into CI band).
#' The last two are each a list with a matrix (return levels)
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2015 + Dec 2016
#' @seealso \code{\link{plotLexBoot}}, \code{\link{distLextreme}}
#' @keywords bootstrap montecarlo hplot dplot distribution ts
#' @export
#' @importFrom berryFunctions logSpaced quantileMean rainbow2
#' 
#' @examples
#' 
#' data(annMax)
#' dlf <- distLextreme(annMax, selection=c("gum","gev","wak","nor"))
#' dlfB <- distLexBoot(dlf, nbest=4, conf.lev=0.5, n=10) # n low for quick example tests
#' plotLexBoot(dlfB)
#' plotLexBoot(dlfB, selection=c("nor","gev"))
#' plotLexBoot(dlfB, selection=c("gum","gev","wak","nor"), order=FALSE)
#' 
#' @param dlf       \code{dlf} object, as returned by \code{\link{distLextreme}}
#' @param nbest     Number of best fitted distribution functions in dlf for which
#'                  bootstrapping is to be done. Overridden by \code{selection}. DEFAULT: 3
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
#' @param quiet     Logical: suppress messages? See \code{\link{distLquantile}}.
#'                  DEFAULT: FALSE
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
quiet=FALSE
)
{
# Selection
if(is.null(selection)) selection <- rownames(dlf$gof)[order(dlf$gof$RMSE)][1:nbest]
# Return period vector:
RPdef <- berryFunctions::logSpaced(min=1, n=100, plot=FALSE,
                                   base=if(log) 1.1708 else 1, max=length(dlf$dat)*2)
if(is.null(RPs)) RPs <-  unique(round(RPdef, digits=2))
# subsample size:
sss <- round(length(dlf$dat_full)*prop)
# Actual computation for each subsample
if(progbars) replicate <- pbapply::pbreplicate
simQ_orig <- replicate(n=n,
   distLquantile(x=c(sample(dlf$dat_full,size=sss)), selection=selection, order=FALSE,
   probs=1-1/(RPs*dlf$npy), empirical=FALSE, progbars=FALSE, time=FALSE, quiet=quiet,
   truncate=dlf$truncate),
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
dlf$exBootCL <- conf.lev
dlf$exBootRPs <- RPs
dlf$distnames <- selection
dlf$distcols <- berryFunctions::rainbow2(length(selection))
dlf$distselector <- "distLexBoot"
invisible(dlf)
}
