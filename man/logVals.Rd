\name{logVals}
\alias{logVals}
\title{Create log-axis values and labels}
\description{create nice values and labels to write at logartihmic axes}
\usage{logVals(from=-7, to=7, Range=range(from, to), base=1,
big.mark="'", decimal.mark=".", scientific=FALSE)}
\arguments{
  \item{from}{lower exponent \emph{OR} vector with data}
  \item{to}{high end}
  \item{Range}{or give from and to as range}
  \item{base}{bases to be used, eg. c(1,2,5)}
  \item{big.mark}{symbol separating thousands, eg. space, comma, dot, etc. see \code{\link{format}} and \code{\link{prettyNum}}}
  \item{decimal.mark}{character separating comma values, see \code{\link{format}} and \code{\link{prettyNum}}}
  \item{scientific}{see \code{\link{format}}}
}
\value{A list with
  \item{vals}{Values for lines and label positions}
  \item{labs}{Formatted values for labels}
  \item{all}{Values for lines}
}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Feb 2014}
\seealso{ \code{\link{log10}}}
\examples{
# Easiest use: vector with data (logVals automatically finds range):
y <- 10^runif(50, -1, 2)
plot(y, log="y") # not much control over placement and format of labels
plot(y, log="y", yaxt="n")
# now do this better, with custom bases:
lv <- logVals(y, base=c(1,2,5) )
axis(2, lv$vals, lv$labs, las=1)

# Default arguments:
lv <- logVals()
str(lv) # values, formatted labels, all 10^x values for lines
plot(1, ylim=c(1e-3, 1e4), log="y", yaxt="n", yaxs="i")
abline(h=lv$all, col=8 )
box("plot")
axis(2, lv$vals, lv$labs, las=1)
lines(seq(0.5, 1.5, len=50), 10^runif(50, -3, 4), col=2)

# Formatting labels:
logVals(                )$labs
logVals(scient=TRUE     )$labs
logVals(big.mark=" "    )$labs
logVals(big=".", dec=",")$labs # German style (not recommended)
}
\keyword{ aplot }
\keyword{ dplot }
