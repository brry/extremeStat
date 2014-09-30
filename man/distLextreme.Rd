\name{distLextreme}
\alias{distLextreme}
\title{Extreme value distributions}
\description{Extreme value statistics for flood risk estimation.
Fits several distributions based on linear moments, calculate return levels for given return periods.
Plots distributions and adds plotting positions by Weibull and Gringorton.
}
\usage{distLextreme(dat, datname, dlf, RPs=c(2,5,10,20,50), speed=TRUE, selection, gofProp=1,
progbars=TRUE, plot=TRUE, add=FALSE, nbest=5, xlim, ylim, las=1,
main=datname, xlab="Return Period RP  [a]",
ylab="Discharge HQ  [m^3/s]", col="black", pch=c(16,3), cex=1,
coldist=rainbow2(nbest), lwd=1, legargs=NULL, ... )}
\arguments{
  \item{dat}{Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.}
  \item{datname}{Character string for main and output list. DEFAULT: internal \code{substitute(dat)}}
  \item{dlf}{List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat!}
  \item{RPs}{ReturnPeriods for which discharge is estimated. DEFAULT: c(2,5,10,20,50)}
  \item{speed}{If TRUE, several computationally intensive time consuming distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}. DEFAULT: TRUE}
  \item{selection}{Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}. Overrides nbest}
  \item{gofProp}{Upper proportion of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail. DEFAULT: 1}
  \item{progbars}{Show progress bars for each loop and \code{\link{message}} execution time? DEFAULT: TRUE}
  \item{plot}{Should the return periods and nbest fitted distributions be plotted? DEFAULT: TRUE}
  \item{add}{If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one. DEFAULT: FALSE}
  \item{nbest}{Number of distributions plotted, in order of goodness of fit. DEFAULT: 5}
  \item{xlim}{X-axis limits. DEFAULT: xlim of plotting positions}
  \item{ylim}{Y-lim. DEFAULT: from min to extended max}
  \item{las}{LabelAxisStyle to orient labels, see \code{\link{par}}. DEFAULT: 1}
  \item{main}{Header of plot. DEFAULT: datname}
  \item{xlab}{X axis label. DEFAULT: "Return Period RP  [a]"}
  \item{ylab}{Y axis label. DEFAULT: "Discharge HQ  [m^3/s]"}
  \item{col}{Plotting point color. DEFAULT: "black"}
  \item{pch}{point characters for plotting positions after Weibull and Grongorton, respectively. DEFAULT: c(16,3)}
  \item{cex}{Character EXpansion of plotting points. DEFAULT: 1}
  \item{coldist}{Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow. DEFAULT: rainbow2(nbest)}
  \item{lwd}{Line WiDth of distribution lines. DEFAULT: 1}
  \item{legargs}{list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd. DEFAULT: NULL}
  \item{\dots}{Further arguments passed to \code{\link{plot}} and \code{\link{lines}}, like lty, lwd, type, pch, ...}
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
# Basic examples
# Goodness of fit
# 1. Measures for GOF
# 2. weights based on GOF
# 3. weighted mean based on GOF
# 4. effect of data proportion used to estimate GOF

# Data from river in Austria
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6, 81.9,
60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1, 40.4,
48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)

# Basic examples
library(lmomco, berryFunctions)
dle <- distLextreme(annMax)
str(dle, max.lev=2)
dle$returnlev

# Estimate discharge that empirically could occur every 80 years:
Q80 <- distLextreme(dlf=dle, RPs=80)$returnlev
sort(Q80[,1], decr=TRUE)
# 99 to 143 m^3/s can make a relevant difference in engineering!

\dontrun{
# compare with other packages:
library2(extRemes)
plot(fevd(annMax))
par(mfrow=c(1,1))
return.level(fevd(annMax, type="GEV")) # "GP", "PP", "Gumbel", "Exponential"
distLextreme(dlf=dle, RPs=c(2,20,100))$returnlev["gev",]
# differences are small, but noticeable...
# if you have time for a more thorough control, please pass me the results!
}

# Line colors / selection distributions to be plotted:
distLextreme(dlf=dle)
distLextreme(dlf=dle, nbest=18, coldist=heat.colors(18))
distLextreme(dlf=dle, selection=c("gev", "gam", "gum"), coldist=3:5)
distLextreme(dlf=dle, legarg=list(cex=0.5, x="bottom", box.col="red", col=3))
# col is ignored. should be specified inside distLextreme
distLextreme(dlf=dle, pch=c(1,NA))

# Goodness of fit
# 1. Measures for GOF ---------------------------------------------
# how well do the distributions fit the original data?
pairs(dle$gof) # measures of goodness of fit are correlated quite well here.
# apparently, it does not really matter which one is used.
dle_nile <- distLextreme(Nile)
pairs(dle_nile$gof)

# compare Ranks by different GOFs
par(mfrow=c(1,2))
distLgofplot(dle)
distLgof(dle, gofProp=0.8, plot=TRUE) # more in examples in ?distLgof
par(mfrow=c(1,1))


# 2. weights based on GOF -----------------------------------------
# normalized Goodness of Fit:
GF <- dle$gof[,"RMSE"] # the lower, the better, the more weight
GF <- max(GF)-GF#+mean(GF) # with mean added, the worst fit is not completely excluded
GF <- GF/sum(GF)
names(GF) <- rownames(dle$gof)
#
GF1 <- dle$gof[,"RMSE"]
GF1 <- max(GF1)-GF1+min(GF1)
GF1 <- GF1/sum(GF1)
#
plot(sort(GF1, decr=TRUE), ylim=c(0, 0.1), type="o", xaxt="n", ylab="Weight",
     las=1, xlab="Dist", main="Weight of distributions for average")
lines(sort(GF, decr=TRUE), pch=3, col=2, type="o")
axis(1, 1:18, names(sort(GF, decr=TRUE)), las=2)
legend("bottomleft", c("max(r)-r + min(r)","max(r)-r"), col=1:2, pch=1:2, lty=1,
       title="Weighted by RMSE (r)")


# 3. weighted mean based on GOF -----------------------------------
distLextreme(dlf=dle, nbest=18, coldist="grey")
# range of discharges:
yval <- seq(from=par("usr")[3], to=par("usr")[4], length=200)
# add weighted average of distributions:
CDFs <- sapply(names(GF), function(d) plmomco(yval,dle$parameter[[d]]))
head(round(CDFs*100), 20) # has NAs out of support region
# remove those completely from weights
GF[ apply(CDFs, 2, function(x) any(is.na(x)))  ]  <- 0
GF <- GF/sum(GF)
weightedMean <- colSums( apply(CDFs, 1, function(x) x*GF), na.rm=TRUE )
lines(1/(1-weightedMean), yval, lwd=3)

# Discharge estimated for 50 years return period, with weighted average:
sum(GF * dle$returnlev[,"RP.50"]) # 115.7
abline(h=115.6561, v=50, col=2)

\dontrun{
# 4. effect of data proportion used to estimate GOF ---------------
Goodness <- function(gofProp)
{
esannMax <- distLextreme(dlf=dle, gofProp=gofProp, plot=FALSE, progbars=FALSE)
esannMax <- cbind(esannMax$returnlev, esannMax$gof)
# normalized Goodness of Fit:
GF <- esannMax[,"RMSE"] # the lower, the better, the more weight
GF <- max(GF)-GF#+mean(GF) # with mean added, the worst fit is not completely excluded
GF <- GF/sum(GF)           # plot(GF1, ylim=c(0, 0.1)); points(GF, pch=3)
# simple mean:              # plot(sort(GF))
av_simple <- mean(esannMax[,"RP.50"])     # 116.8383
# Discharge estimated for 50 years return period, bei weighted average:
av_weight <- sum(GF * esannMax[,"RP.50"]) # 117.2784
# mean of best 3 extreme value distribution functions:
av_3best <- mean(esannMax[ order(esannMax[,"RMSE"])[1:3] , "RP.50"]) # 118.3722
# most functions underestimate the discharge, if we assume that the weibull
# PP method correctly calculates the return period of the highest value
c(av_simple=av_simple, av_weight=av_weight, av_3best=av_3best)
}

Goodness(0.1)
Goodness(0.2)
Proportion <- seq(0.05, 1, len=100)
GoF <- sapply(Proportion, Goodness)
plot(Proportion, GoF[3,], type="l")
lines(Proportion, GoF[2,], col=2)
lines(Proportion, GoF[1,], col=3)
title(main="The proportion of the data included for calculating RMSE does not matter!")
# but it did before, so there may be an error here
}

}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
\keyword{ts}

