\name{distLextremePlot}
\alias{distLextremePlot}
\title{Plot extreme value statistics}
\description{Plots distributions fitted by L-moments and adds plotting positions by Weibull and Gringorton.}
\usage{distLextremePlot(dlf, selection, add=FALSE, nbest=5, xlim=NULL, ylim=NULL,
las=1, main=dlf$datname, xlab="Return Period RP  [a]", ylab="Discharge HQ  [m^3/s]",
col="black", pch=c(16,3), cex=1, coldist=rainbow2(nbest), lwd=1,
legend=TRUE, legargs=NULL, linargs=NULL, quiet=FALSE, ... )}
\arguments{
  \item{dlf}{List as returned by \code{\link{distLextreme}}, containing the elements \code{dat, parameter, gof}. Overrides dat!}
  \item{selection}{Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}}
  \item{add}{If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one. DEFAULT: FALSE}
  \item{nbest}{Number of distributions plotted, in order of goodness of fit. DEFAULT: 5}
  \item{xlim}{X-axis limits. DEFAULT: xlim of plotting positions}
  \item{ylim}{Y-lim. DEFAULT: from min to extended max}
  \item{las}{LabelAxisStyle to orient labels, see \code{\link{par}}. DEFAULT: 1}
  \item{main}{Title of plot. DEFAULT: dlf$datname}
  \item{xlab}{X axis label. DEFAULT: "Return Period RP  [a]"}
  \item{ylab}{Y axis label. DEFAULT: "Discharge HQ  [m^3/s]"}
  \item{col}{Plotting point color. DEFAULT: "black"}
  \item{pch}{point characters for plotting positions after Weibull and Grongorton, respectively. NA to suppress in plot and legend. DEFAULT: c(16,3)}
  \item{cex}{Character EXpansion of plotting points. DEFAULT: 1}
  \item{coldist}{Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow}
  \item{lwd}{Line WiDth of distribution lines. DEFAULT: 1}
  \item{legend}{Logical. Add a legend? DEFAULT: TRUE}
  \item{legargs}{list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd. DEFAULT: NULL}
  \item{linargs}{list of arguments passed to \code{\link{lines}} like lty, lwd, type, pch, ... DEFAULT: NULL}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link{plot}}, like log="x", xaxt="n", ...}
}
\details{This is an auxiliary graphing function to \code{\link{distLextreme}}}
\value{none, plots things.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, March 2015}
\seealso{\code{\link{distLextreme}}, \code{\link{distLplot}}  }
\examples{
#see
?distLextreme
}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
