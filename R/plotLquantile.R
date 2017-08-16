#' Plot quantiles of distributions fitted with L-moments
#' 
#' @return invisible dlf object, see \code{\link{printL}}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLquantile}}, \code{\link{plotLfit}}
#' @keywords hplot distribution
#' @export
#' @importFrom berryFunctions rainbow2 owa
#' @importFrom graphics lines par
#' @importFrom stats quantile
#' @examples
#' # See distLquantile
#' 
#' @param dlf     List as returned by \code{\link{distLquantile}}, containing the
#'                elements \code{dat, parameter, gof, datname, quant}
#' @param nbest,selection,order Distributions to be plotted, see \code{\link{plotLfit}}
#' @param rows    Rowname(s) of \code{dlf$quant} that should be drawn instead of
#'                the selection / nbest highest ranking distribution functions.
#'                'GPD*' will select all the gpd fits. heights and distcols must then
#'                accordingly have at least 13 elements (or will be recycled).
#'                DEFAULT: NULL
#' @param heights Coordinates of quantile line ends, recycled if necessary.
#'                DEFAULT: 20\% of plot height.
#' @param distcols Color for each distribution added with \code{\link{lines}}.
#'                DEFAULT: dlfplot$distcols
#' @param linargs Arguments passed to \code{\link{lines}}. DEFAULT: NULL
#' @param \dots   Further arguments passed to \code{\link{plotLfit}}
#' 
plotLquantile <- function(
dlf,
nbest=5,
selection=NULL,
order=FALSE,
rows=NULL,
heights=stats::quantile(par("usr")[3:4], 0.2),
distcols=dlfplot$distcols,
linargs=NULL,
... )
{
if(is.null(dlf$quant)) stop("dlf does not contain 'quant' element.")
dlfplot <- plotLfit(dlf, nbest=nbest, selection=selection, order=order, ...)
dlfplot$distselector <- "plotLquantile"

# For which distributions should quantile lines be plotted:
dn <- dlfplot$distnames
if(!is.null(rows))
  {
  rnq <- rownames(dlf$quant)
  if(any(rows=="GPD*")) rows <- c(rows[rows!="GPD*"], rnq[grepl("GPD_",rnq)])
  dn <- unique(rows)
  dlfplot$distnames <- dn
  dlfplot$distcols <- berryFunctions::rainbow2(length(dn))
  }

# prepare lines
heights <- rep(heights, len=length(dn))
distcols <- rep(distcols, len=length(dn)) # recycle 2
columns <- colnames(dlf$quant)
columns <- columns[grepl("%", columns)]

# actually add lines
for(i in columns) do.call(graphics::lines, args=berryFunctions::owa(list(
      x=dlf$quant[dn,i], y=heights, type="h", col=distcols), linargs, "x","y"))

invisible(dlfplot)
}
