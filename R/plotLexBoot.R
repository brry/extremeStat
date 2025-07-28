#' Bootstrapping uncertainty intervals for return periods
#' 
#' plot bootstrap uncertainty intervals for \code{\link{plotLextreme}}.
#' 
#' @return invisible dlf object, see \code{\link{printL}}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLexBoot}}
#' @keywords bootstrap montecarlo hplot dplot distribution ts
#' @export
#' @importFrom berryFunctions ciBand
#' 
#' @examples
#' # see distLexBoot
#' 
#' @param dlf       \code{dlf} object, as returned by \code{\link{distLexBoot}}
#' @param selection Character vector with distribution function names to be used.
#'                  Suggested to keep this low. DEFAULT: NULL
#' @param add       Add to existing plot? DEFAULT: FALSE
#' @param log       Plot on a logarithmic axis. DEFAULT: TRUE
#' @param \dots     Further arguments passed to \code{\link{plotLextreme}}.
#'                  If add=TRUE, they are instead passed to
#'                  \code{berryFunctions::\link[berryFunctions]{ciBand}}
#' 
plotLexBoot <- function(
dlf,
selection=NULL,
add=FALSE,
log=TRUE,
...
)
{
# Selection
if(is.null(selection)) selection <- dlf$distnames
# plot
if(!add) dlf <- plotLextreme(dlf=dlf, selection=selection, log=log, ...)
exBootCI <- dlf$exBootCI[selection]
for(i in length(exBootCI):1) if(add)
berryFunctions::ciBand(yu=exBootCI[[i]][2,], yl=exBootCI[[i]][1,], x=dlf$exBootRPs,
                       add=TRUE, colm=dlf$distcols[i], nastars=FALSE, ...) else
berryFunctions::ciBand(yu=exBootCI[[i]][2,], yl=exBootCI[[i]][1,], x=dlf$exBootRPs,
                       add=TRUE, colm=dlf$distcols[i], nastars=FALSE)
# Output
invisible(dlf)
}
