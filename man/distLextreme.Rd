\name{distLextreme}
\alias{distLextreme}
\title{Extreme value distributions}
\description{Extreme value statistics for flood risk estimation.
Fits several distributions based on linear moments, calculate return levels for given return periods.
Plots distributions and adds plotting positions by Weibull and Gringorton.
}
\usage{distLextreme(dat, datname, dlf, RPs=c(2,5,10,20,50), speed=TRUE, ks=FALSE,
selection, order=FALSE, gofProp=1, progbars=length(dat)>200, time=TRUE,
plot=TRUE, add=FALSE,
nbest=5, xlim, ylim, las=1, main=datname, xlab="Return Period RP  [a]",
ylab="Discharge HQ  [m\\U00B3/s]", col="black", pch=c(16,3), cex=1,
coldist=rainbow2(nbest), lwd=1, legend=TRUE, legargs=NULL, linargs=NULL,
quiet=FALSE, ... )}
\arguments{
  \item{dat}{Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.}
  \item{datname}{Character string for main and output list. DEFAULT: internal \code{substitute(dat)}}
  \item{dlf}{List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat!}
  \item{RPs}{ReturnPeriods for which discharge is estimated. DEFAULT: c(2,5,10,20,50)}
  \item{speed}{If TRUE, several computationally intensive time consuming distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}. DEFAULT: TRUE}
  \item{ks}{Include ks.test results in \code{dlf$gof}? Computing is much faster on big datasets when FALSE. DEFAULT: TRUE}
  \item{selection}{Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}. Overrides nbest}
  \item{order}{If selection is given, should legend and colors in \code{\link{distLextremePlot}} be ordered by gof anyways?}
\item{gofProp}{Upper proportion of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail. DEFAULT: 1}
  \item{progbars}{Show progress bars for each loop and \code{\link{message}} execution time? DEFAULT: TRUE if n > 200}
  \item{time}{\code{\link{message}} execution time? DEFAULT: TRUE}
  \item{plot}{Should the return periods and nbest fitted distributions be plotted? DEFAULT: TRUE}
  \item{add}{If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one. DEFAULT: FALSE}
  \item{nbest}{Number of distributions plotted, in order of goodness of fit. DEFAULT: 5}
  \item{xlim}{X-axis limits. DEFAULT: xlim of plotting positions}
  \item{ylim}{Y-lim. DEFAULT: from min to extended max}
  \item{las}{LabelAxisStyle to orient labels, see \code{\link{par}}. DEFAULT: 1}
  \item{main}{Header of plot. DEFAULT: datname}
  \item{xlab}{X axis label. DEFAULT: "Return Period RP  [a]"}
  \item{ylab}{Y axis label. Please note that the ubuntu pdf viewer might be unable to display unicode superscript. DEFAULT: "Discharge HQ  [m\\U00B3/s]"}
  \item{col}{Plotting point color. DEFAULT: "black"}
  \item{pch}{point characters for plotting positions after Weibull and Grongorton, respectively. NA to suppress in plot and legend. DEFAULT: c(16,3)}
  \item{cex}{Character EXpansion of plotting points. DEFAULT: 1}
  \item{coldist}{Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow. DEFAULT: rainbow2(nbest)}
  \item{lwd}{Line WiDth of distribution lines. DEFAULT: 1}
  \item{legend}{Logical. Add a legend? DEFAULT: TRUE}
  \item{legargs}{list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd. DEFAULT: NULL}
  \item{linargs}{list of arguments passed to \code{\link{lines}} like lty, lwd, type, pch, ... DEFAULT: NULL}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link{plot}}, like log="x", xaxt="n", ...}
}
\details{
this adds weibull and gringorton plotting positions to the distribution lines, which are estimated from the linear moments of the data itself.\cr
I personally believe that if you have, say, 35 values in \code{dat}, the highest return period should be around 36 years (Weibull) and not 60 (Gringorton).\cr
The plotting positions don't affect the distribution parameter estimation, so this dispute is not really important.
But if you care, go ahead and google "weibull vs gringorton plotting positions".
}
\value{List as explained in \code{\link{extremeStat}}. The added element is \code{returnlev}, a data.frame with the return level (discharge) for all given RPs and for each distribution.}
\note{This function replaces \code{berryFunctions::extremeStatLmom}}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, 2012-2014}
\references{\url{http://RclickHandbuch.wordpress.com} Chapter 15 (German)\cr
Christoph Mudersbach: Untersuchungen zur Ermittlung von hydrologischen Bemessungsgroessen mit Verfahren der instationaeren Extremwertstatistik }
\seealso{\code{\link{distLfit}}. \code{\link[extRemes]{fevd}} in the package \code{extRemes}.}
\examples{
# Data from river in Austria
# Basic examples
# Advanced options
# weighted mean based on Goodness of fit (GOF)
# Effect of data proportion used to estimate GOF
# compare extremeStat with other packages

library(lmomco)
library(berryFunctions)

# Data from river in Austria
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6, 81.9,
60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1, 40.4,
48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)

# Basic examples ---------------------------------------------------------------
dle <- distLextreme(annMax, log="x")
abline(h=1:10*20, col=8) ;  berryFunctions::logAxis(1)
dle <- distLextreme(annMax, nbest=10, add=TRUE, legargs=list(bg="white", cex=0.9))
# Why do I not get the following warning in interactive mode?
## calculating RMSE:
## Warning in rmse(tcdfs[[d]], ecdfs, quiet = q
##   7 NAs were omitted from 35 data points.

# Object structure:
str(dle, max.lev=2)
distLprint(dle)

# discharge levels for default return periods:
dle$returnlev

# Estimate discharge that could occur every 80 years (at least empirically):
Q80 <- distLextreme(dlf=dle, RPs=80)$returnlev
round(sort(Q80[,1], decr=TRUE),1)
# 99 to 143 m^3/s can make a relevant difference in engineering!

# Return period of a given discharge value:
sort(1/(1-sapply(dle$parameter, plmomco, x=150) ))
# exponential: 105 years, gev (general extreme value dist): 536,
# Weibull: every 1305 years only
# averages of distribution functions weighted by fit quality: see below



# Advanced options -------------------------------------------------------------
distLextremePlot(dlf=dle)
# Line colors / select distributions to be plotted:
distLextremePlot(dle, nbest=17, coldist=heat.colors(17), lty=1:5) # lty is recycled
distLextremePlot(dle, selection=c("gev", "gam", "gum"), coldist=3:5, col=3, lty=3:2)
# Why do I not get the following warning in interactive mode?
## Warning in if (is.na(dn)) stop("No distributions are left with selection ",  :
## the condition has length > 1 and only the first element will be used
# apparently, warnings do not get passed from one function to the next...

distLextremePlot(dle, legarg=list(cex=0.5, x="bottom", box.col="red", col=3))
# col in legarg list is (correctly) ignored
distLextremePlot(dle, pch=c(1,NA)) # only Weibull plotting positions
# add different dataset to existing plot:
distLextreme(Nile/15, add=TRUE, pch=NA, coldist=1, selection="wak", legend=FALSE)

# Logarithmic axis
distLextreme(Nile, log="x", xaxt="n", pch=NA, coldist=NA, legend=FALSE)
logAxis(1)
distLextreme(Nile, add=TRUE, legargs=c(bg="white"), nbest=8)



# Weighted mean based on GOF ---------------------------------------------------
# GOF: Goodness of fit, for details see
?distLgof
# weights:
distLgofPlot(dle, ranks=FALSE)

# Discharge estimated for 50 years return period, with weighted average:
sum(dle$gof$weight1 * dle$returnlev[,"RP.50"]) # 114.9
sum(dle$gof$weight2 * dle$returnlev[,"RP.50"]) # 115.3
sum(dle$gof$weight3 * dle$returnlev[,"RP.50"]) # 115.9
distLextremePlot(dle, nbest=17, legend=FALSE)
abline(h=115.9, v=50)

# discharge estimate continuously across all RPs:
yval <- seq(from=par("usr")[3], to=par("usr")[4], length=200)  # range of discharges
# add weighted average of distributions:
CDFs <- sapply(rownames(dle$gof), function(d) plmomco(yval,dle$parameter[[d]]))
weightedMean <- colSums( apply(CDFs, 1, function(x) x * dle$gof$weight3), na.rm=TRUE )
lines(1/(1-weightedMean), yval, lwd=3)


head(round(CDFs*100), 20) # gov has NAs out of support region
# remove those completely from weights
dle$gof$weight1[ apply(CDFs, 2, function(x) any(is.na(x)))  ]  <- 0
dle$gof$weight1 <- dle$gof$weight1/sum(dle$gof$weight1)
weightedMean <- colSums( apply(CDFs, 1, function(x) x * dle$gof$weight1), na.rm=TRUE )
lines(1/(1-weightedMean), yval, lwd=3, col="orange") # not very different



# Effect of data proportion used to estimate GOF -------------------------------
\dontrun{ # removed from package testing due to computing time
# Discharge estimated for 50 years return period
Goodness <- function(gofProp)
{ # ExtremeStatistics
ES <- distLextreme(annMax, gofProp=gofProp, plot=FALSE, progbars=FALSE)
ES <- cbind(ES$returnlev, ES$gof)
# simple mean:              # plot(sort(GF))
av_simple <- mean(ES[,"RP.50"])     # 114.0178  old: 116.8383
# weighted average:
av_weight <- sum(ES$weight2 * ES[,"RP.50"]) # 115.2842  old: 117.2784
# mean of best 3 distribution functions:
av_3best <- mean(ES[1:3, "RP.50"]) # 114.7579   old: 118.3722
# most functions underestimate the discharge, if we assume that the weibull
# PP method correctly calculates the return period of the highest value
c(av_simple=av_simple, av_weight=av_weight, av_3best=av_3best)
}

Goodness(1)
Goodness(0.2)
Proportion <- seq(0.05, 1, len=50)
GoF <- sapply(Proportion, Goodness) # takes a while to compute
plot(Proportion, GoF[3,], type="l", las=1)
lines(Proportion, GoF[2,], col=2)
lines(Proportion, GoF[1,], col=3)
legend("bottomright", col=1:3, lty=1,
   legend=c("mean of 3 dists with lowest RMSE", "weighted average", "simple average"))
title(main="The proportion of the data included\nfor calculating RMSE does matter!")



# compare extremeStat with other packages: ---------------------------------------
library(extRemes)
plot(fevd(annMax))
par(mfrow=c(1,1))
return.level(fevd(annMax, type="GEV")) # "GP", "PP", "Gumbel", "Exponential"
distLextreme(dlf=dle, RPs=c(2,20,100))$returnlev["gev",]
# differences are small, but noticeable...
# if you have time for a more thorough control, please pass me the results!
}

# yet another dataset for testing purposes:
Dresden_AnnualMax <- c(403, 468, 497, 539, 542, 634, 662, 765, 834, 847, 851, 873,
885, 983, 996, 1020, 1028, 1090, 1096, 1110, 1173, 1180, 1180,
1220, 1270, 1285, 1329, 1360, 1360, 1387, 1401, 1410, 1410, 1456,
1556, 1580, 1610, 1630, 1680, 1734, 1740, 1748, 1780, 1800, 1820,
1896, 1962, 2000, 2010, 2238, 2270, 2860, 4500)
distLextreme(Dresden_AnnualMax)

}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
\keyword{ts}
