#' Plot distributions fitted with linear moments
#' 
#' Plot histogram and distribution densities \emph{or} ecdf with cumulated probability
#' 
#' @details By default, this plots density instead of CDF, because the distributions are
#' easier to discern and tail behaviour is easier to judge visually. See also
#' urlhttp://www.vosesoftware.com/vosesoftware/ModelRiskHelp/index.htm#Presenting_results/Cumulative_plots/Relationship_between_cdf_and_density_\%28histogram\%29_plots.htm
#' 
#' @return dlf with coldist added, returned invisibly.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{distLfit}}, \code{\link{distLextreme}}
#' @keywords hplot distribution
#' @export
#' @importFrom berryFunctions lim0 owa rainbow2
#' @importFrom lmomco plmomco dlmomco supdist
#' @examples
#'  # See distLfit and distLquantile
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof, datname, gofProp}
#' @param nbest Number of distributions plotted, in order of goodness of fit. DEFAULT: 5
#' @param selection Names of distributions in \code{dlf$parameter} that will be drawn. Overrides nbest. DEFAULT: NULL
#' @param order If selection is given, should legend and colors be ordered by gof anyways? DEFAULT: FALSE
#' @param cdf If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE
#' @param log If TRUE, logAxis is called. DEFAULT: FALSE
#' @param percentline If TRUE, draw vertical line at 1-dlf$gofProp of dlf$dat. If NA, only do so if gofProp!=1. DEFAULT: NA
#' @param percentargs List of arguments passed to \code{\link{abline}}. DEFAULT: NULL
#' @param supportends If TRUE, dots are placed at the support bounds. DEFAULT: TRUE
#' @param breaks \code{\link{hist}} breaks. DEFAULT: 20
#' @param xlim \code{\link{hist}} or \code{\link{ecdf}} xlim. DEFAULT: extendrange(dat, f=0.15)
#' @param ylim \code{\link{hist}} or \code{\link{ecdf}} ylim. DEFAULT: NULL
#' @param xaxs,yaxs \code{\link{hist}} or \code{\link{ecdf}} xaxs and yaxs. DEFAULT: both "i"
#' @param xaxt \code{\link{par}} xaxt. "n" suppresses axis and numbers, which is used if log
#' @param col \code{\link{hist}} bar color or \code{\link{ecdf}} point color. DEFAULT: "grey"
#' @param main,xlab,ylab \code{\link{hist}} or \code{\link{ecdf}} main, xlab, ylab. DEFAULT: internal abstraction fom \code{dlf$datname}
#' @param las Label Axis Style for orientation of numbers along axes. DEFAULT: 1
#' @param coldist Color for each distribution added with \code{\link{lines}}. DEFAULT: \code{\link[berryFunctions]{rainbow2}}
#' @param lty Line TYpe for plotted distributions. Recycled vector of length nbest. DEFAULT: 1
#' @param add If TRUE, hist is not called before adding lines. This lets you add lines highly customized one by one. DEFAULT: FALSE
#' @param logargs List of arguments passed to \code{\link{logAxis}} if \code{log=TRUE}. DEFAULT: NULL
#' @param legend Should \code{\link{legend}} be called? DEFAULT: TRUE
#' @param legargs List of arguments passed to \code{\link{legend}} except for legend and col. DEFAULT: NULL
#' @param histargs List of arguments passed to \code{\link{hist}} or \code{\link{ecdf}} except for x, freq. DEFAULT: NULL
#' @param qlines Should vertical lines marking the quantiles be added?
#'        Ignored if \code{dlf} does not contain the element \code{quant}. DEFAULT: FALSE
#' @param qheights Coordinates of quantile line ends, recycled if necessary. DEFAULT: 20\% of plot height.
#' @param qrow Rowname(s) of \code{dlf$quant} that should be drawn instead of the nbest highest ranking distribution functions.
#'             'q_gpd*' will select all the gpd fits, including gpa.
#'             qheights and coldist must then accordingly have at least 13 elements (or will be recycled). DEFAULT: NULL
#' @param qlinargs Arguments passed to \code{\link{lines}} for qlines. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{lines}}, like type, pch, ...
#' 
distLplot <- function(
dlf,
nbest=5,
selection=NULL,
order=FALSE,
cdf=FALSE,
log=FALSE,
percentline=NA,
percentargs=NULL,
supportends=TRUE,
breaks=20,
xlim=extendrange(dlf$dat, f=0.15),
ylim=NULL,
xaxs="i", yaxs="i",
xaxt,
col="grey",
main, xlab, ylab,
las=1,
coldist=berryFunctions::rainbow2(nbest),
lty=1,
add=FALSE,
logargs=NULL,
legend=TRUE,
legargs=NULL,
histargs=NULL,
qlines=FALSE,
qheights=stats::quantile(par("usr")[3:4], 0.2),
qrow=NULL,
qlinargs=NULL,
... )
{
# input checks:
if(!is.list(dlf)) stop("dlf must be a list.")
# checking list elements:
if(is.null(dlf$dat)) stop("dlf must contain the element dat")
if(is.null(dlf$parameter)) stop("dlf must contain the element parameter")
if(is.null(dlf$gof)) stop("dlf must contain the element gof")
if(is.null(dlf$datname)) stop("dlf must contain the element datname")
if(is.null(dlf$gofProp)) stop("dlf must contain the element gofProp")
# distribution selection:
if(!is.null(selection))
  {
  names(selection) <- selection
  if(order) selection <- selection[rownames(dlf$gof)]
  selection <- selection[!is.na(selection)]
  sing <- selection %in% rownames(dlf$gof)
  if(!any(sing)) stop("selection (", pastec(selection[!sing]), ") is not available in dlf$gof.")
  if(any(!sing)) 
    {
    curselsing <- pastec(selection[!sing])
    on.exit(message("Note in disLplot: selection (", curselsing, ") is not available in dlf$gof."), add=TRUE)
    }
  selection <- selection[sing]
  dlf$gof <- dlf$gof[selection, ]      # continuing with a different variant of gof!!
  nbest <- length(selection)
  }
# input checks:
if(nbest < 0) stop("nbest must be positive")
if(nbest > nrow(dlf$gof) & is.null(qrow)) {nbest <- nrow(dlf$gof)}
# internal defaults:
if(missing(xaxt)       ) xaxt <- if(log) "n" else "s"
if(missing(xlab)       ) xlab <- dlf$datname
#
if(!is.null(qrow))
  {
  if(any(qrow=="q_gpd*")) qrow <- c(qrow[qrow!="q_gpd*"], "gpa", "q_gpd_evir_pwm",
         "q_gpd_evir_ml", "q_gpd_evd", "q_gpd_extRemes_MLE", "q_gpd_extRemes_GMLE",
         "q_gpd_extRemes_Bayesian", "q_gpd_extRemes_Lmoments", "q_gpd_fExtremes_pwm",
         "q_gpd_fExtremes_mle", "q_gpd_ismev", "q_gpd_Renext_r", "q_gpd_Renext_f")
  dnqrow <- unique(qrow)
  }
#
# draw histogram or ecdf -------------------------------------------------------
if(!add)
{
  if(cdf)
  {
  if(is.null(ylim)) ylim <- c(dlf$truncate,1)
  if(missing(ylab)) ylab <- "(Empirical) Cumulated Density (CDF)"
  if(missing(main)) main <- paste("Cumulated density distributions of", dlf$datname)
  ecdfdef <- list(x=ecdf(dlf$dat_full), do.points=TRUE, col=col, xlim=xlim, xaxt=xaxt, ylab=ylab,
             ylim=ylim, xaxs=xaxs, yaxs=yaxs, main=main, xlab=xlab, las=las)
  do.call(plot, args=berryFunctions::owa(ecdfdef, histargs, "x", "y"))
  if(log)
    {do.call(logAxis, args=berryFunctions::owa(list(xaxt="s"), logargs))
     do.call(lines,   args=berryFunctions::owa(ecdfdef,       histargs, "x", "y"))
    }
  }
  else # if not cdf, then density
  {
  if(is.null(ylim)) ylim <- berryFunctions::lim0(hist(dlf$dat, breaks=breaks, plot=FALSE)$density,
                                  curtail=if(yaxs=="i") FALSE else TRUE)
  if(missing(ylab)) ylab <- "Probability Density Function (PDF)"
  if(missing(main)) main <- paste("Density distributions of", dlf$datname)
  op <- par(xaxs=xaxs, yaxs=yaxs, xaxt=xaxt)
  histdef <- list(x=dlf$dat, breaks=breaks, col=col, xlim=xlim, ylim=ylim, ylab=ylab,
                  freq=FALSE, main=main, xlab=xlab, las=las)
  do.call(hist, args=berryFunctions::owa(histdef, histargs, "x", "freq"))
  if(log)
    {do.call(logAxis, args=berryFunctions::owa(list(xaxt="s"), logargs))
     do.call(hist,    args=berryFunctions::owa(c(histdef, add=TRUE), histargs, "x", "freq", "add"))
    }
  par(op)
  }
}
#
# distribution names and colors: -----------------------------------------------
if(nbest < 1) return(invisible(NULL)) # and stop executing
dn <- rownames(dlf$gof)[1:nbest]
if(qlines) if(!is.null(qrow)) dn <- dnqrow[dnqrow %in% rownames(dlf$gof)]
coldist <- rep(coldist, length=nbest)  # recycle 1
if(length(dlf$parameter)<1)
  {
  message("Note in distLplot: no distributions were available in dlf$parameter.")
  dlf$coldist <- coldist
  return(invisible(dlf))
  }
lty <- rep(lty, length=nbest)
# add distribution function lines:
if(cdf) lfun <- lmomco::plmomco else lfun <- lmomco::dlmomco
for(i in length(dn):1)
  {
  xval <- seq(from=par("usr")[1], to=par("usr")[2], length=300)
  # cut xval to support region
  paramd <- dlf$parameter[[dn[i]]]
  #browser()
  if(is.null(paramd)) next # 
  if(!is.null(paramd$ifail)) if(paramd$ifail!=0) next # emu # this may be too strict...
  if(paramd$type=="gld") next # lmomco 2.2.2 cdfgld bug
  if(!is.null(paramd$error)) if(paramd$convergence!=0) next # gep 
  xsup <- lmomco::supdist(paramd)$support
  xval <- xval[ xval>xsup[1] & xval<xsup[2] ]
  # only plot distribution line if there is some support:
  if(length(xval)>0){
  yval <- lfun(xval,paramd)
  if(cdf & dlf$truncate!=0) yval <- yval*(1-dlf$truncate) + dlf$truncate ## yval <- (yval-dlf$truncate)/(1-dlf$truncate)
  lines(xval, yval, col=coldist[i], lty=lty[i], ...)
  if(supportends)
    {
    # last point within support range, if support ends in graphing region:
    lo <- if(xsup[1] > par("usr")[1])      xval[1] else NA
    hi <- if(xsup[2] < par("usr")[2]) tail(xval,1) else NA
    if(!is.na(lo) ) points(lo,      yval[1], pch=16, col=coldist[i])
    if(!is.na(hi) ) points(hi, tail(yval,1), pch=16, col=coldist[i])
    } # end if supportends
  } # end if xval has values
  } # end for loop over distribution functions
# draw vertical gofProp line:
if(is.na(percentline)) percentline <- if(dlf$gofProp!=1) TRUE else FALSE
if(percentline) do.call(abline, args=berryFunctions::owa(list(
       v=quantile(dlf$dat, probs=1-dlf$gofProp), lty=3, col="red"), percentargs))
# draw quantile lines:
if(qlines)
  {
  if(is.null(dlf$quant)) stop("qlines=TRUE, but dlf does not contain 'quant' element.")
  if(!is.null(qrow)) dn <- dnqrow
  # qheights <- lfun(x=qval, para=dlf$parameter[[use[i]]])
  qheights <- rep(qheights, len=length(dn))
  coldist <- rep(coldist, len=length(dn)) # recycle 2
  for(i in seq_len(ncol(dlf$quant))) do.call(graphics::lines, args=berryFunctions::owa(list(
      x=dlf$quant[dn,i], y=qheights, type="h", col=coldist), qlinargs, "x","y"))
  }
# legend - write the names of distributions:
legdef <- list(legend=dn, lwd=1, col=coldist, x="right", cex=0.7, lty=lty)
if(legend) do.call(graphics::legend, args=berryFunctions::owa(legdef, legargs, "legend","col","lty"))
# add to (or change) output:
dlf$coldist <- coldist
return(invisible(dlf))
} # end function
