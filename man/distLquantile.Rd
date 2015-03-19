\name{distLquantile}
\alias{distLquantile}
\title{distribution quantiles}
\description{Parametrical quantiles of distributions fitted to a sample.}
\usage{distLquantile(x=NULL, probs=0:4/4, truncate=0, type=NULL, dlf=NULL, order=TRUE,
plot=FALSE, cdf=FALSE, lines=TRUE, linargs=NULL, empirical=TRUE, ... )}
\arguments{
  \item{x}{Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used. DEFAULT: NULL}
  \item{probs}{Numeric vector of probabilities with values in [0,1]. DEFAULT: 0:4/4}
  \item{truncate}{Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly. DEFAULT: 0}
  \item{type}{Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used. DEFAULT: NULL}
  \item{dlf}{dlf object described in \code{\link{extremeStat}}. Ignored if x is given. DEFAULT: NULL}
  \item{order}{Sort results by GOF? If TRUE (the default) and length(type)>1, the output is ordered by dlf$gof, else by order of appearance in dlf$parameter. DEFAULT: TRUE}
  \item{plot}{Should \code{\link{distLplot}} be called? DEFAULT: FALSE}
  \item{cdf}{If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE}
  \item{lines}{Should vertical lines marking the quantiles be added? DEFAULT: TRUE}
  \item{linargs}{Arguments passed to \code{\link{lines}}. DEFAULT: NULL}
  \item{empirical}{Add vertical line for empirical \code{\link{quantileMean}}? DEFAULT: TRUE}
  \item{\dots}{Arguments passed to \code{\link{distLfit}} if x or truncate is given. Passed to \code{\link{distLplot}} if plot=TRUE and x is not given and truncate is 0.}
}
\details{
Very high quantiles (99\% and higher) need large sample sizes for \code{\link{quantile}} to yield a robust estimate.
Theoretically, at least 1/(1-probs) values must be present, e.g. 10'000 for Q99.99\%.
With smaller sample sizes (eg n=35), they underestimate the actual (but unknown) quantile.
Parametrical quantiles need only small sample sizes. They don't have a systematical underestimation bias, but have higher variability.
}
\value{Matrix with distribution quantile values. Returns NAs for probs below truncate.}
\note{NAs are always removed from x in \code{\link{distLfit}}}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, March 2015}
\references{My master thesis paper, to be submitted soon.}
\seealso{\code{\link{distLfit}},
Xian Zhou, Liuquan Sun and Haobo Ren (2000): Quantile estimation for left truncated and right censored data, Statistica Sinica 10
\url{http://www3.stat.sinica.edu.tw/statistica/oldpdf/A10n411.pdf}\cr
require("truncdist") }
\examples{
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6,
81.9, 60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1,
40.4, 48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)
dlf <- distLfit(annMax)

distLquantile(annMax) # fits parameters. If speed is important in big datasets
# and parameters are already available, pass them via dlf:
distLquantile(dlf=dlf, probs=0:10/10, type=c("wak", "gev", "kap"), order=FALSE)

distLquantile(annMax, plot=FALSE, type="wak")

# censored (truncated, trimmed) quantile:
qwak <- distLquantile(annMax, type="wak") ; abline(v=qwak, col=2, lwd=3)
qwak2 <- distLquantile(annMax, type="wak", truncate=0.2, add=TRUE, coldist=4)
abline(v=qwak2, col=4, lwd=2)

# Simulation of truncation effect
library(lmomco)
#set.seed(42)
rnum <- rlmomco(n=1e3, para=dlf$parameter$gev)
myprobs <- c(0.9, 0.95, 0.99, 0.999)
mytrunc <- seq(0, 0.9, length.out=20)
trunceffect <- sapply(mytrunc, function(mt) distLquantile(rnum, type="gev",
                      probs=myprobs, truncate=mt, plot=FALSE, progbars=FALSE)  )
#
op <- par(mfrow=c(2,1), mar=c(2,4.5,2,0.5), cex.main=1)
distLquantile(rnum, type="gev", progbars=FALSE, ylab="", xlab="")
distLquantile(rnum, type="gev", truncate=0.3, add=TRUE, progbars=FALSE, coldist=4)
legend("right", c("fitted GEV", "fitted with truncate=0.3"), lty=1, col=c(2,4), bg="white")
par(mar=c(3,4.5,3,0.5))
plot(mytrunc, trunceffect[1,], ylim=range(trunceffect), las=1, type="l",
     main=c("High quantiles of 1000 random numbers from gev distribution",
           "Estimation based on proportion of lower values truncated"),
     xlab="", ylab="parametrical quantile")
title(xlab="Proportion censored", mgp=c(1.8,1,0))
for(i in 2:4) lines(mytrunc, trunceffect[i,])
library("berryFunctions")
textField(rep(0.5,4), trunceffect[,11], paste0("Q",myprobs*100,"\%") )
par(op)
}
\keyword{distribution}
\keyword{robust}
\keyword{univar}
