# distributions via linear moments
# Berry Boessenkool, Sept 2014
# Plot density distributions fitted by distLfit

distLplot <- function(
dlf,                    # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof, datname, gofProp}
nbest=5,                # Number of distributions plotted, in order of goodness of fit
selection=NULL,         # Names of distributions in \code{dlf$parameter} that will be drawn. Overrides nbest.
order=FALSE,            # If selection is given, should legend and colors be ordered by gof anyways?
cdf=FALSE,              # If TRUE, plot cumulated DF instead of probability density
log=FALSE,              # If TRUE, logAxis is called.
percentline=NA,         # If TRUE, draw vertical line at 1-dlf$gofProp of dlf$dat. If NA, only do so if gofProp!=1
percentargs=NULL,       # List of arguments passed to \code{\link{abline}}.
supportends=TRUE,       # If TRUE, dots are placed at the support bounds.
breaks=20,              # \code{\link{hist}} breaks
xlim=extendrange(dat, f=0.15), # \code{\link{hist}} or \code{\link{ecdf}} xlim
ylim=NULL,              # \code{\link{hist}} or \code{\link{ecdf}} ylim
xaxs="i", yaxs="i",     # \code{\link{hist}} or \code{\link{ecdf}} xaxs and yaxs. DEFAULT: both "i"
xaxt,                   # \code{\link{par}} xaxt. "n" suppresses axis and numbers, which is used if log
col="grey",             # \code{\link{hist}} bar color or \code{\link{ecdf}} point color
main, xlab, ylab,       # \code{\link{hist}} or \code{\link{ecdf}} main, xlab, ylab
las=1,                  # Label Axis Style for orientation of numbers along axes.
coldist=rainbow2(nbest),# Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow2
add=FALSE,              # If TRUE, hist is not called before adding lines. This lets you add lines highly customized one by one
logargs=NULL,           # List of arguments passed to \code{\link{logAxis}} if \code{log=TRUE}
legend=TRUE,            # Should \code{\link{legend}} be called?
legargs=NULL,           # List of arguments passed to \code{\link{legend}} except for legend and col
histargs=NULL,          # List of arguments passed to \code{\link{hist}} or \code{\link{ecdf}} except for x, freq
... )                   # Further arguments passed to \code{\link{lines}}, like lty, type, pch, ...
{
# input checks:
if(!is.list(dlf)) stop("dlf must be a list.")
# Objects from list:
dat       <- dlf$dat
parameter <- dlf$parameter
gof       <- dlf$gof
# checking list elements:
if(is.null(dat)) stop("dlf must contain the element dat")
if(is.null(parameter)) stop("dlf must contain the element parameter")
if(is.null(gof)) stop("dlf must contain the element gof")
if(is.null(dlf$datname)) stop("dlf must contain the element datname")
if(is.null(dlf$gofProp)) stop("dlf must contain the element gofProp")
# distribution selection:
if(!is.null(selection))
  {
  names(selection) <- selection
  if(order) selection <- selection[rownames(gof)]
  selection <- selection[!is.na(selection)]
  sing <- selection %in% rownames(gof)
  if(!any(sing)) stop("selection ", paste(selection[!sing], collapse=", "), " is not available in dlf$gof.")
  selection <- selection[sing]
  gof <- gof[selection, ]
  nbest <- length(selection)
  }
# input checks:
if(nbest < 0) stop("nbest must be positive")
if(nbest > nrow(gof)) {nbest <- nrow(gof)}
# internal defaults:
if(missing(xaxt)       ) xaxt <- if(log) "n" else "s"
if(missing(xlab)       ) xlab <- dlf$datname
# draw histogram or ecdf
if(!add)
  if(cdf)
  {
  if(is.null(ylim)) ylim <- c(0,1)
  if(missing(ylab)) ylab <- "(Empirical) Cumulated Density (CDF)"
  if(missing(main)) main <- paste("Cumulated density distributions of", dlf$datname)
  ecdfdef <- list(x=ecdf(dat), col=col, xlim=xlim, xaxt=xaxt, ylab=ylab,
             ylim=ylim, xaxs=xaxs, yaxs=yaxs, main=main, xlab=xlab, las=las)
  do.call(plot, args=owa(ecdfdef, histargs, "x", "y"))
  if(log)
    {do.call(logAxis, args=owa(list(xaxt="s"), logargs))
     do.call(lines,   args=owa(ecdfdef,       histargs, "x", "y"))
    }
  }
  else # if not cdf, then density
  {
  if(is.null(ylim)) ylim <- lim0(hist(dat, breaks=breaks, plot=FALSE)$density,
                                  curtail=if(yaxs=="i") FALSE else TRUE)
  if(missing(ylab)) ylab <- "Probability Density Function (PDF)"
  if(missing(main)) main <- paste("Density distributions of", dlf$datname)
  op <- par(xaxs=xaxs, yaxs=yaxs, xaxt=xaxt)
  histdef <- list(x=dat, breaks=breaks, col=col, xlim=xlim, ylim=ylim, ylab=ylab,
                  freq=FALSE, main=main, xlab=xlab, las=las)
  do.call(hist, args=owa(histdef, histargs, "x", "freq"))
  if(log)
    {do.call(logAxis, args=owa(list(xaxt="s"), logargs))
     do.call(hist,    args=owa(c(histdef, add=TRUE), histargs, "x", "freq", "add"))
    }
  par(op)
  }
# distribution names and colors:
if(nbest < 1) return(invisible(NULL)) # and stop executing
dn <- rownames(gof)[1:nbest]
coldist <- rep(coldist, length=nbest)
# add distributions:
if(cdf) lfun <- plmomco else lfun <- dlmomco
for(i in nbest:1)
  {
  xval <- seq(from=par("usr")[1], to=par("usr")[2], length=300)
  # cut xval to support region
  paramd <- parameter[[dn[i]]]
  xsup <- supdist(paramd)$support
  xval <- xval[ xval>xsup[1] & xval<xsup[2] ]
  # last point within support range, if support ends in graphing region:
  lo <- if(xsup[1] > par("usr")[1])      xval[1] else NA
  hi <- if(xsup[2] < par("usr")[2]) tail(xval,1) else NA
  # only plot distribution line if there is some support:
  if(length(xval)>0)            lines(xval, lfun(xval,paramd), col=coldist[i], ...)
  if(supportends & !is.na(lo) ) points(lo, lfun(lo,paramd), pch=16, col=coldist[i])
  if(supportends & !is.na(hi) ) points(hi, lfun(hi,paramd), pch=16, col=coldist[i])
  }
# draw vertical gofProp line:
if(is.na(percentline)) percentline <- if(dlf$gofProp!=1) TRUE else FALSE
if(percentline) do.call(abline, args=owa(list(v=quantile(dat, probs=1-dlf$gofProp),
                                                lty=3, col="red"), percentargs))
# legend - write the names of distributions:
legdef <- list(legend=dn, lwd=1, col=coldist, x="right", cex=0.7)
if(legend) do.call(graphics::legend, args=owa(legdef, legargs, c("legend","col")))
# add to output:
dlf$coldist <- coldist
return(invisible(dlf))
} # end function
