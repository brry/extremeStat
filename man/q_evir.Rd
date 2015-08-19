\name{q_evir}
\alias{q_evir}
\alias{q_evir2}
\title{GPD quantile of sample}
\description{Compute quantile of General Pareto Distribution fitted to sample
by peak over treshold (POT) method declaring treshold from truncation proportion
using the established evir (extreme values in r) package}
\usage{q_evir(x, probs, truncate, pngdev=TRUE, quiet=FALSE, ...)}
\arguments{
  \item{x}{Vector with values}
  \item{probs}{Probabilities of truncated (Peak over treshold) quantile}
  \item{truncate}{Truncation percentage (proportion of sample discarded)}
  \item{pngdev}{sink \code{evir::quant} graph output to file (is removed later) instead of openening \code{\link{dev.new}}, which also is closed later. Using TRUE avoids the graphics device showing, but \code{getwd()} must be writable. DEFAULT: TRUE}
  \item{quiet}{Should messages be suppressed? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link[evir]{quant}}}
}
\details{Opens and closes external device to suppress \code{evir::quant} plotting.\cr
q_evir2 doesn't do that, but is a manual and slightly changed copy from the computing part of the code}
\value{Vector of quantile estimates for each value of probs}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, July 2015}
\references{\url{http://stats.stackexchange.com/questions/129438/different-quantiles-of-a-fitted-gpd-in-different-r-packages}}
\seealso{\code{\link{distLquantile}}, \code{evir::\link[evir]{quant}} }
\examples{
library(lmomco)
library(berryFunctions)
library(evir)

# I. GPD is not necessarily applicable for full samples!
set.seed(42); rnum <- rbeta(1000, 2, 7)*300
d <- distLquantile(rnum, plot=TRUE, probs=0.99, nbest=18   )
# gpa (=GPD) performs badly here, compared to the other distributions, because:
distLgofPlot(distLfit(rnum, quiet=TRUE), ranks=FALSE) # GPA rmse is relatively high...
d <- distLquantile(rnum, plot=TRUE, probs=0.99, selection=c("wei", "gpa"),
              xlim=c(120, 270), ylim=c(0, 0.003), quiet=TRUE)
# and GPA has limited distribution support!! (dot at end)
# (computing quantiles or drawing random numbers will never exceed 184)
# Remember: the dist. parameters are in distLfit(rnum, sel="gpa")$param$gpa
par <- lmomco::pargpa(lmomco::lmoms(rnum))
lmomco::quagpa(f=1, para=par)


# II. GPA performs better for distribution tails: 2nd extreme value theorem
browseURL("https://en.wikipedia.org/wiki/Pickands-Balkema-de_Haan_theorem")
distLquantile(rnum, plot=TRUE, probs=0.99, progbars=FALSE)
distLquantile(rnum, plot=TRUE, probs=0.99, truncate=0.8, progbars=FALSE)
distLquantile(rnum, plot=TRUE, probs=0.99, selection=c("wak", "gpa"),
              xlim=c(120, 270), ylim=c(0.6, 1), truncate=0.8, cdf=TRUE)


evir::quant(rnum, p=0.95, models=70, end=1000)
abline(v=-c(100,200), lty=2) # stable quantile estimate at the tail, but not for very high truncation
abline(h=q_evir(rnum, probs=0.95, truncate=0.8), lty=2)


# III. evir::quant vs linear moments lmomco (implemented in extremeStat::distLquantile)
# a) nontruncated full sample:
comparequantiles <- function(trunc=0, xlim=c(0,210), ylim=c(trunc,1), ... )
{
probs <- seq(trunc,1, length=200)
d <- distLquantile(rnum, probs=probs, truncate=trunc, ...)
plot(ecdf(rnum), las=1, ylim=ylim, col=8, col.01line=NA, xlim=xlim, 
     xlab="Quantile estimate", main=paste(trunc*100, "% truncated"))
lines(d["quantileMean",], probs, col="orange")
lines(d["weighted1",], probs, col=1)
lines(d["gpa",],       probs, col=2)
lines(d["q_evir",],    probs, col=4)
legend("bottomright", c("ecdf","quantileMean","weighted","gpa","evir"), 
       col=c(8,"orange",1,2,4), lty=1)
invisible(d)
}
comparequantiles(ylim=lim0(1), xlim=c(0,210) )
abline(h=c(0.3,0.9), lty=2)
# evir GPD method underestimates 30% Quantile, overestimates at 90%
# lmomco gpa does the same, but not so extremely.
# Hardly surprising: GPD is not a beta distibution!

# Here are the numbers
dlq <- distLquantile(rnum, probs=c(0.3, 0.9), quiet=TRUE)
quantileMean(rnum, probs=c(0.3, 0.9))       # empirical
dlq["weighted1",]                           # fitted distributions
dlq["gpa",]                                 # GPD fitted by linear moments
q_evir(rnum, probs=c(0.3, 0.9), truncate=0) # GPD in package evir

         
# Now let's look at the tail of the beta-distributed random numbers:
# b) truncated Peak Over Treshold (POT) EV
comparequantiles(trunc=0.8, ylim=c(0.95, 1), xlim=c(140,210) )
abline(h=0.99, lty=2) 
# lmomco and evir yield the same result. Both slightly overestimate

comparequantiles(trunc=0.9, ylim=c(0.95, 1), xlim=c(140,210), quiet=TRUE )
d <- comparequantiles(trunc=0.98, xlim=c(155,255), quiet=TRUE )
}
\keyword{distribution}
\keyword{robust}
\keyword{univar}
