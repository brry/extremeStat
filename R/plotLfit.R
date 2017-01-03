#' Plot distributions fitted with L-moments
#' 
#' Plot histogram and distribution densities \emph{or} ecdf with cumulated probability
#' 
#' @details By default, this plots density instead of CDF, because the distributions are
#' easier to discern and tail behaviour is easier to judge visually. See also
#' \url{http://www.vosesoftware.com/ModelRiskHelp/index.htm#Presenting_results/Cumulative_plots/Relationship_between_cdf_and_density_(histogram)_plots.htm}
#' 
#' @return invisible dlf object, see \code{\link{printL}}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{distLfit}}, \code{\link{plotLquantile}}
#' @keywords hplot distribution
#' @export
#' @importFrom berryFunctions lim0 owa rainbow2
#' @importFrom lmomco plmomco dlmomco supdist
#' @examples
#'  # See distLfit
#' 
#' @param dlf        List as returned by \code{\link{distLfit}}, containing the 
#'                   elements \code{dat, parameter, gof, datname}
#' @param nbest      Number of distributions plotted, in order of goodness of fit. 
#'                   DEFAULT: 5
#' @param selection  Names of distributions in \code{dlf$parameter} that will be drawn. 
#'                   Overrides nbest. DEFAULT: NULL
#' @param order      Logical: order legend and colors by RMSE, even if dlf$gof is 
#'                   unordered or selection is given? DEFAULT: TRUE
#' @param rmse       Integers. If rmse != 0, RMSE values are added to legend.
#'                   They are rounded to \code{rmse} digits. DEFAULT: 4
#' @param cdf        If TRUE, plot cumulated DF instead of probability density. 
#'                   DEFAULT: FALSE
#' @param log        If TRUE, logAxis is called. DEFAULT: FALSE
#' @param supportends If TRUE, dots are placed at the support bounds. DEFAULT: TRUE
#' @param breaks     \code{\link{hist}} breaks. DEFAULT: 20
#' @param xlim,ylim  \code{\link{hist}} or \code{\link{ecdf}} axis limits.
#' @param xaxs,yaxs  \code{\link{hist}} or \code{\link{ecdf}} xaxs and yaxs. DEFAULT: both "i"
#' @param xaxt       \code{\link{par}} xaxt. "n" suppresses axis and numbers, 
#'                   which is used if log
#' @param col        \code{\link{hist}} bar color or \code{\link{ecdf}} point color. 
#'                   DEFAULT: "grey"
#' @param main,xlab,ylab \code{\link{hist}} or \code{\link{ecdf}} main, xlab, ylab. 
#'                   DEFAULT: abstractions fom \code{dlf$datname}
#' @param las        Label Axis Style for orientation of numbers along axes. DEFAULT: 1
#' @param distcols   Color for each distribution added with \code{\link{lines}}. 
#'                   DEFAULT: \code{\link[berryFunctions]{rainbow2}}
#' @param lty        Line TYpe for plotted distributions. 
#'                   Recycled vector of length nbest. DEFAULT: 1
#' @param add        If TRUE, hist/ecdf is not called before adding lines. 
#'                   This lets you add lines highly customized one by one. 
#'                   DEFAULT: FALSE
#' @param logargs    List of arguments passed to \code{\link{logAxis}} if 
#'                   \code{log=TRUE}. DEFAULT: NULL
#' @param legend     Should \code{\link{legend}} be called? DEFAULT: TRUE
#' @param legargs    List of arguments passed to \code{\link{legend}} except for 
#'                   legend and col. DEFAULT: NULL
#' @param histargs   List of arguments passed to \code{\link{hist}} or \code{\link{ecdf}} 
#'                   except for x, freq. DEFAULT: NULL
#' @param \dots      Further arguments passed to \code{\link{lines}}, like type, pch, ...
#' 
plotLfit <- function(
dlf,
nbest=5,
selection=NULL,
order=TRUE,
rmse=4,
cdf=FALSE,
log=FALSE,
supportends=TRUE,
breaks=20,
xlim=extendrange(dlf$dat, f=0.15),
ylim=NULL,
xaxs="i", yaxs="i",
xaxt=if(log) "n" else "s",
col="grey",
main=paste(if(cdf)"Cumulated", "density distributions of", dlf$datname), 
xlab=dlf$datname, 
ylab=if(cdf) "(Empirical) Cumulated Density (CDF)" else "Probability Density Function (PDF)",
las=1,
distcols=berryFunctions::rainbow2(nbest),
lty=1,
add=FALSE,
logargs=NULL,
legend=TRUE,
legargs=NULL,
histargs=NULL,
... )
{
# input checks:
if(!is.list(dlf)) stop("dlf must be a list.")
# checking list elements:
if(is.null(dlf$dat))       stop("dlf must contain the element dat")
if(is.null(dlf$parameter)) stop("dlf must contain the element parameter")
if(is.null(dlf$gof))       stop("dlf must contain the element gof")
if(is.null(dlf$datname))   stop("dlf must contain the element datname")
# distribution selection:
dn <- rownames(dlf$gof)[!is.na(dlf$gof$RMSE)]
if(!is.null(selection))
  {
  names(selection) <- selection
  sing <- selection %in% dn
  if(!any(sing)) stop("selection (", toString(selection[!sing]), ") is not available in dlf.")
  if(any(!sing)) message("Note in disLplot: selection (",toString(selection[!sing]),
                         ") is not available in dlf.")
  selection <- selection[sing]
  nbest <- length(selection)
  dn <- selection
  }
# final selection + nbest checks:
nfitted <- length(dn)
if(nbest > nfitted) {nbest <- nfitted}
if(nbest < 1) stop("nbest must be a positive integer > 0")
if(order) dn <- dn[order(dlf$gof[dn,"RMSE"])]
dn <- dn[1:nbest]
distcols <- rep(distcols, length=nbest)  # recycle
lty <- rep(lty, length=nbest)
#
#
# draw histogram or ecdf -------------------------------------------------------
if(!add)
{
  if(cdf)
  {
  if(is.null(ylim)) ylim <- c(dlf$truncate,1)
  ecdfdef <- list(x=ecdf(dlf$dat_full), do.points=TRUE, col=col, xlim=xlim, xaxt=xaxt, 
             ylab=ylab, ylim=ylim, xaxs=xaxs, yaxs=yaxs, main=main, xlab=xlab, las=las)
  do.call(plot, args=berryFunctions::owa(ecdfdef, histargs, "x", "y"))
  if(log)
    {do.call(logAxis, args=berryFunctions::owa(list(xaxt="s"), logargs))
     do.call(lines,   args=berryFunctions::owa(ecdfdef,       histargs, "x", "y"))
    }
  }
  else # if not cdf, then density
  {
  if(is.null(ylim)) ylim <- berryFunctions::lim0(hist(dlf$dat, breaks=breaks,plot=FALSE)$density,
                                  curtail=if(yaxs=="i") FALSE else TRUE)
  op <- par(xaxs=xaxs, yaxs=yaxs, xaxt=xaxt)
  histdef <- list(x=dlf$dat, breaks=breaks, col=col, xlim=xlim, ylim=ylim, ylab=ylab,
                  freq=FALSE, main=main, xlab=xlab, las=las)
  do.call(hist, args=berryFunctions::owa(histdef, histargs, "x", "freq"))
  if(log)
    {do.call(logAxis, args=berryFunctions::owa(list(xaxt="s"), logargs))
     do.call(hist,    args=berryFunctions::owa(c(histdef, add=TRUE), histargs, 
                                               "x", "freq", "add"))
    }
  par(op)
  }
}
#
# add distribution function lines: ---------------------------------------------
if(cdf) lfun <- lmomco::plmomco else lfun <- lmomco::dlmomco
for(i in length(dn):1)
  {
  xval <- seq(from=par("usr")[1], to=par("usr")[2], length=300)
  # cut xval to support region
  paramd <- dlf$parameter[[dn[i]]]
  xsup <- lmomco::supdist(paramd)$support
  xval <- xval[ xval>xsup[1] & xval<xsup[2] ]
  # only plot distribution line if there is some support:
  if(length(xval)>0){
  yval <- lfun(xval,paramd)
  if(cdf & dlf$truncate!=0) yval <- yval*(1-dlf$truncate) + dlf$truncate 
                         ## yval <- (yval-dlf$truncate)/(1-dlf$truncate)
  lines(xval, yval, col=distcols[i], lty=lty[i], ...)
  if(supportends)
    {
    # last point within support range, if support ends in graphing region:
    lo <- if(xsup[1] > par("usr")[1])      xval[1] else NA
    hi <- if(xsup[2] < par("usr")[2]) tail(xval,1) else NA
    if(!is.na(lo) ) points(lo,      yval[1], pch=16, col=distcols[i])
    if(!is.na(hi) ) points(hi, tail(yval,1), pch=16, col=distcols[i])
    } # end if supportends
  } # end if xval has values
  } # end for loop over distribution functions
#
# legend -----------------------------------------------------------------------
legdn <- paste(if(rmse[1]>0) formatC(round(dlf$gof[dn,"RMSE"],rmse), format='f', digits=rmse), dn)
legdef <- list(legend=legdn, lwd=1, col=distcols, x="right", cex=0.7, lty=lty,
               title="distribution GOF")
if(legend) do.call(graphics::legend, args=berryFunctions::owa(legdef, legargs, 
                                                              "legend","col","lty"))
# add to (or change) output:
dlf$distnames <- dn
dlf$distcols <- distcols
dlf$distselector <- "plotLfit"
return(invisible(dlf))
} # end function
