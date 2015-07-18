\name{extremeStat}
\alias{extremeStat-package}
\alias{extremeStat}
\docType{package}
\title{Extreme value statistics on a linear scale}
\description{
fit (via linear moments), plot (on a linear scale) and compare (by goodness of fit)
several (extreme value) distributions to estimate discharge at given return periods.\cr
This package heavily relies on and thankfully acknowledges the package lmomco by WH Asquith.
}
\details{The common object to share between functions is a list (\code{dlf}) with:\cr
\tabular{ll}{
\code{dat}       \tab numeric vector with (extreme) values \cr
\code{datname}   \tab character string for main, xlab etc \cr
\code{gofProp}   \tab number between 0 and 1; upper proportion of \code{dat} to compute goodness of fit from\cr
\code{parameter} \tab list (usually of length 17 if \code{speed=TRUE}) with parameters of each distribution\cr
\code{gof}       \tab dataframe with 'Goodness of Fit' measures, sorted by RMSE of theoretical and empirical cumulated density\cr
\code{returnlev} \tab dataframe with values of distributions for given return periods (\code{RPs}). This element is only added in \code{\link{distLextreme}}\cr
\code{RP___}     \tab Return periods according to plotting positions, see below. \cr
\code{coldist}   \tab Colors for plotting, added in \code{\link{distLplot}}\cr
\code{truncate}  \tab Truncation percentage, only relevant for \code{\link{distLquantile}}\cr
}
It can be printed with \code{\link{distLprint}}, which may be transformed to a class with printing method.\cr
PP:\cr
Plotting positions are not used for fitting distributions, but for plotting only\cr
The ranks of ascendingly sorted extreme values are used to compute the probability of non-exceedence Pn:\cr
\code{Pn_w <-  Rank      /(n+1)       # Weibull}
\code{Pn_g <- (Rank-0.44)/(n+0.12)    # Gringorton (taken from lmom:::evplot.default)}\cr
Finally: RP = Returnperiod = recurrence interval = 1/P_exceedence = 1/(1-P_nonexc.), thus:\cr
\code{RPweibull = 1/(1-Pn_w)} and analogous for gringorton.\cr

The main functions in the extremeStat package are:
\tabular{ll}{
\code{\link{distLextreme}}     \tab analyse extreme value statistics, calls \code{distLfit} and \code{distLextremePlot}.\cr
\code{\link{distLextremePlot}} \tab plot distribution lines and plotting positions.\cr
\code{\link{distLfit}}         \tab fit the parameters, calls \code{gof} and \code{distLplot}.\cr
\code{\link{distLgof}}         \tab calculate goodness of fits, calls \code{distLgofPlot}. Can also be executed with \code{dlf} to minimize computing time by not fitting the parameters again.\cr
\code{\link{distLplot}}        \tab plot density or cumulated density of data and distributions.\cr
\code{\link{distLgofPlot}}     \tab compare distribution ranks of different \code{distLgof} methods.\cr
}
\code{Depends} on 'berryFunctions' for \code{\link{rmse}}, \code{\link{rsquare}}, \code{\link{logAxis}}, \code{\link{logVals}}.\cr
\code{Suggests} 'pbapply' to see progress bars if you have large (n > 1e3) datasets.\cr
At some places you'll find \code{## not run} in the examples.
These code blocks were excluded from checking while building,
mainly because they are computationally intensive and should not take so much of CRANs resources.
Normally, you should be able to run them in an interactive session.\cr
If you do find unexecutable code, please tell me!\cr
This package was motivated by my need to compare the fits of several distributions to data.
It was originally triggered by a flood estimation assignment we had in class 2012,
and it bothered me that we just assumed the gumbel distribution would fit the data fine.\cr
With the updated form of the original function, I think this is a useful package to compare fits.\cr
I'm no expert on distributions, so I welcome all suggestions you might have for me.
}
\seealso{
If you are looking for more detailed (uncertainty) analysis, eg confidence intervals,
check out the package \code{extRemes}, especially the function \code{\link[extRemes]{fevd}}.
\url{http://cran.r-project.org/package=extRemes}\cr
Intro slides: \url{http://sites.lsa.umich.edu/eva2015/wp-content/uploads/sites/44/2015/06/Intro2EVT.pdf}\cr
Parameter fitting and distribution functions: \url{http://cran.r-project.org/package=lmomco}\cr
Distributions: \url{https://www.rmetrics.org/files/Meielisalp2009/Presentations/Scott.pdf}
and: \url{http://cran.r-project.org/web/views/Distributions.html} \cr
R in Hydrology: \url{http://abouthydrology.blogspot.de/2012/08/r-resources-for-hydrologists.html}\cr
}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, 2014-2015}
\examples{
# annual discharge maxima from a stream in Austria
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6,
   81.9, 60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1,
   40.4, 48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)
# time plot:
plot(annMax, type="l")
require("berryFunctions"); lines(movAv(annMax), lwd=2, col=2)
dle <- distLextreme(annMax)
cbind(dle$returnlev, dle$gof)
}
\keyword{package}
\keyword{documentation}
