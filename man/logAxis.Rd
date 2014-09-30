\name{logAxis}
\alias{logAxis}
\title{label logarithmic axes}
\description{shortcut to calling \code{\link{logVals}}, \code{\link{axis}} and \code{\link{abline}}}
\usage{logAxis(side=1, log=NULL, lcol="grey", lty=1, lwd=1, expr, las=1, from=-7,
to=7, Range=range(from, to), base=c(1,2,5), big.mark="'", decimal.mark=".",
scientific=FALSE, ...)}
\arguments{
  \item{side}{Which \code{\link{axis}} is to be labeled?. DEFAULT: 1}
  \item{log}{Is the axis logarithmic by plot(log="x")? internal DEFAULT: par("xlog") or "ylog". DEFAULT: NULL}
  \item{lcol}{Color of gridlines drawn in the graph with \code{\link{abline}}, NA to suppress. DEFAULT: "grey"}
  \item{lty, lwd}{Type of gridlines. DEFAULT: 1}
  \item{expr}{Expression drawing over the ablines, like (points(x,y). Can be code within {braces}.}
  \item{las}{LabelAxisStyle for the orientation of the labels. DEFAULT: 1}
  \item{from}{Lower exponent OR vector with data, as in \code{\link{logVals}}. DEFAULT: -7}
  \item{to}{High end exponent. DEFAULT: 7}
  \item{Range}{Override from and to as range. DEFAULT: range(from, to)}
  \item{base}{Bases to be used, eg. c(1,2,5) or 1. DEFAULT: c(1,2,5)}
  \item{big.mark}{Symbol separating thousands, eg. space, comma, dot, etc. see "format" and "prettyNum". DEFAULT: "'"}
  \item{decimal.mark}{Character separating comma values, see "format" and "prettyNum". DEFAULT: "."}
  \item{scientific}{See \code{\link{format}}. DEFAULT: FALSE}
  \item{\dots}{further arguments passed to axis, like \code{lwd, col.ticks, hadj, lty}, ...}
}
\value{An invisible list with
  \item{vals}{Values for lines and label positions}
  \item{labs}{Formatted values for labels}
  \item{all}{Values for lines}
}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{logVals}}, \code{\link{log10}}}
\examples{
x <- 10^runif(200, -1, 2)
plot(x, 1:200, xaxt="n", log="x")
logAxis(base=1)
# overplot vertical lines:
logAxis(expr=points(x, 1:200), base=c(1,2,5), col=2 )

# plots where log="x" is not possible:
hist(log10(x), breaks=20, col.axis="grey", main="")
logAxis(side=3, expr=hist(log10(x),breaks=20,add=TRUE, col=3))
}
\keyword{aplot}
\keyword{dplot}
