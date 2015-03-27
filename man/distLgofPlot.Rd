\name{distLgofPlot}
\alias{distLgofPlot}
\title{rank comparison}
\description{plot rank comparison of fitted distributions calculated by \code{\link{distLgof}}.}
\usage{distLgofPlot(dlf, nbest=NA, order=TRUE, add=FALSE,
type="o", col=c(1,3,4), pch=1:3, lty=1, quiet=FALSE, ...)}
\arguments{
  \item{dlf}{List as returned by \code{\link{distLfit}}, containing the elements \code{gof, gofProp}}
  \item{nbest}{Number of distributions plotted by \code{\link{distLplot}}. DEFAULT: NA}
  \item{order}{If TRUE, the best distribution is at the top. False to compare across different dlfs. DEFAULT: TRUE}
  \item{add}{If TRUE, just the lines and legend are added to an existing plot. DEFAULT: FALSE}
  \item{type}{type of line-point combination as in \code{\link{plot}}. DEFAULT: "o"}
  \item{col, pch, lty}{Vectors with 3 values for line customization. Recycled if necessary. DEFAULTS: see created plot.}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link{plot}}.}
}
\value{None.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{distLgof}}, \code{\link{distLfit}}}
\examples{# see distLgof
}
\keyword{hplot}
\keyword{distribution}
