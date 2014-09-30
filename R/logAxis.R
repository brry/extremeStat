# Function to label logarithmic axes
# berry-b@gmx.de, Feb 2014, idea 2013

logAxis <- function(
  side=1,      # Which \code{\link{axis}} is to be labeled?
  log=NULL,    # Is the axis logarithmic by plot(log="x")? internal DEFAULT: par("xlog") or "ylog"
  lcol="grey", # Color of gridlines drawn in the graph with \code{\link{abline}}, NA to suppress.
  lty=1, lwd=1,# Type of gridlines
  expr,        # Expression drawing over the ablines, like (points(x,y). Can be code within {braces}.
  las=1,       # LabelAxisStyle for the orientation of the labels
  from=-7,               # Lower exponent OR vector with data, as in \code{\link{logVals}}
  to=7,                  # High end exponent
  Range=range(from, to), # Override from and to as range
  base=c(1,2,5),         # Bases to be used, eg. c(1,2,5) or 1
  big.mark="'",          # Symbol separating thousands, eg. space, comma, dot, etc. see "format" and "prettyNum"
  decimal.mark=".",      # Character separating comma values, see "format" and "prettyNum"
  scientific=FALSE,      # See \code{\link{format}}
  ...)         # further arguments passed to axis, like \code{lwd, col.ticks, hadj, lty}, ...
{
# get labels and positions:
lv <- logVals(from=from, to=to, Range=Range, base=base, big.mark=big.mark,
              decimal.mark=decimal.mark, scientific=scientific)
# choose vertical or horizontal lines:
if(side==1 | side==3)
  {
  if(is.null(log)) log <- par("xlog")
  if(log) abline(v=lv$all,     col=lcol, lty=lty, lwd=lwd)
  else abline(v=log10(lv$all), col=lcol, lty=lty, lwd=lwd)
  }
else
  {
  if(is.null(log)) log <- par("ylog")
  if(log) abline(h=lv$all,     col=lcol, lty=lty, lwd=lwd)
  else abline(h=log10(lv$all), col=lcol, lty=lty, lwd=lwd)
  }
box("plot")
# axis labels:
if(log) axis(side=side, at=lv$vals,        labels=lv$labs, las=las, ...)
else    axis(side=side, at=log10(lv$vals), labels=lv$labs, las=las, ...)
# overplot ablines with expr:
if(!missing(expr)) expr
# output:
return(invisible(lv))
}
