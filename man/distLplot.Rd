\name{distLplot}
\alias{distLplot}
\title{plot distributions fitted with linear moments}
\description{plot histogram and distribution densities \emph{or} ecdf with cumulated probability}
\usage{distLplot(dlf, nbest=5, selection=NULL, cdf=FALSE, log=FALSE,
percentline=NA, percentargs=NULL, supportends=TRUE,
breaks=20, xlim=extendrange(dat, f=0.15), ylim=NULL, xaxs="i", yaxs="i", xaxt, col="grey",
main, xlab, ylab, las=1, coldist=rainbow2(nbest), add=FALSE,
logargs=NULL, legend=TRUE, legargs=NULL, histargs=NULL, quiet=FALSE, ... )}
\arguments{
  \item{dlf}{List as returned by \code{\link{distLfit}}, containing the elements \code{dat, datname, parameter, gof}}
  \item{nbest}{Number of distributions plotted, in order of goodness of fit. DEFAULT: 5}
  \item{selection}{Names of distributions in \code{dlf$parameter} that will be drawn. Overrides nbest.}
  \item{cdf}{If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE}
  \item{log}{If TRUE, logAxis is called. DEFAULT: FALSE}
  \item{percentline}{If TRUE, draw vertical line at 1-dlf$gofProp of dlf$dat. If NA, only do so if gofProp!=1. DEFAULT: NA}
  \item{percentargs}{List of arguments passed to \code{\link{abline}}  DEFAULT: NULL}
  \item{supportends}{If TRUE, dots are placed at the support bounds. DEFAULT: TRUE}
  \item{breaks}{\code{\link{hist}} breaks. DEFAULT: 20}
  \item{xlim}{\code{\link{hist}} or \code{\link{ecdf}} xlim. DEFAULT: extendrange(dat, f=0.15)}
  \item{ylim}{\code{\link{hist}} or \code{\link{ecdf}} ylim. DEFAULT: NULL}
  \item{xaxs, yaxs}{\code{\link{hist}} or \code{\link{ecdf}} xaxs and yaxs. DEFAULT: both "i"}
  \item{xaxt}{\code{\link{par}} xaxt. "n" suppresses axis and numbers, which is used if log.}
  \item{col}{\code{\link{hist}} bar color or \code{\link{ecdf}} point color. DEFAULT: "grey"}
  \item{main, xlab, ylab}{\code{\link{hist}} or \code{\link{ecdf}} main, xlab, ylab. DEFAULT: internal abstraction fom \code{dlf$datname}}
  \item{las}{Label Axis Style for orientation of numbers along axes. DEFAULT: 1}
  \item{coldist}{Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow2(nbest)}
  \item{add}{If TRUE, hist is not called before adding lines. This lets you add lines highly customized one by one. DEFAULT: FALSE}
  \item{logargs}{List of arguments passed to \code{\link{logAxis}} if \code{log=TRUE}. DEFAULT: NULL}
  \item{legend}{Should \code{\link{legend}} be called? DEFAULT: TRUE}
  \item{legargs}{List of arguments passed to \code{\link{legend}} except for legend and col. DEFAULT: NULL}
  \item{histargs}{List of arguments passed to \code{\link{hist}} or \code{\link{ecdf}} except for x, freq. DEFAULT: NULL}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link{lines}}, like lty, type, pch, ...}
}
\value{dlf with coldist added, returned invisibly.}
\details{By default, this plots density instead of CDF, because the distributions
are easier to discern and tail behaviour is easier to judge visually.
See also url{http://www.vosesoftware.com/vosesoftware/ModelRiskHelp/index.htm#Presenting_results/Cumulative_plots/Relationship_between_cdf_and_density_\%28histogram\%29_plots.htm}
}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{distLfit}}, \code{\link{distLextreme}}}
\examples{ # See distLfit
}
\keyword{hplot}
\keyword{distribution}
