\name{distLgof}
\alias{distLgof}
\title{quality of distribution fits}
\description{calculate goodness of fit for several distributions, plot rank comparison.}
\usage{distLgof(dlf, gofProp, plot=TRUE, progbars=TRUE)}
\arguments{
  \item{dlf}{List as returned by \code{\link{distLfit}}, containing the elements \code{dat, datname, gofProp, parameter}}
  \item{gofProp}{Overrides value in list. Proportion (0:1) of highest values in \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail}
  \item{plot}{Call \code{\link{distLgofplot}}? DEFAULT: TRUE}
  \item{progbars}{Show progress bars for each loop? DEFAULT: TRUE)}
}
\value{List as explained in \code{\link{extremeStat}}. The added element is gof, a data.frame the root mean square error (RMSE) and R squared (R2) for the top \code{gofProp} of \code{dat}, as well as the p and D values from a simple \code{\link{ks.test}} for each distribution function.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{distLgofplot}}, \code{\link{distLfit}}. More complex estimates of quality of fits: \url{http://chjs.soche.cl/papers/vol4n1_2013/ChJS-04-01-04.pdf} }
\examples{
require(lmomco)
exdat <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6, 81.9,
60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1, 40.4,
48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)

# Kolmogorov-Smirnov Tests return slightly different values:
ks.test(exdat, "pnorm", mean(exdat), sd(exdat) )$p.value
ks.test(exdat, "cdfnor", parnor(lmoms(exdat)))$p.value

# effect of Proportion of values used to calculate RMSE
dlf100 <- distLfit(exdat, gofProp=1) # the default gofProp: 100%
distLplot(dlf100, nbest=19, breaks=10, coldist=grey(1:19/19) )
distLgofplot(dlf100)

dlf50 <- distLgof(dlf100, gofProp=0.5)
# so revgum, nor and rice do well on the upper half by R2, but bad by RMSE
distLplot(dlf50, breaks=10)
# The red dashed line shows the cut above which the data were used to get R2/rmse

distLplot(dlf=dlf50, cdf=TRUE)
distLplot(dlf=dlf50, selection=c("pe3","wei", "rice", "nor", "revgum"),
          xlim=c(60,120), ylim=c(0.5, 1), cdf=TRUE, col=1)
dlf50$gof

library(berryFunctions) # for linReg
compranks <- function(d)
{
gofProp <- 0.5
dat2 <- sort(exdat, decreasing=TRUE)[  1:(gofProp*length(exdat))  ]
tcdfs <- plmomco(dat2,dlf50$parameter[[d]])
ecdfs <- ecdf(exdat)(dat2) # Empirical CDF
# Root Mean Square Error, R squared:
linReg(tcdfs, ecdfs, lwd=1, pch=16, main=d, xlim=c(0.5, 1), ylim=c(0.5, 1), pos1="topleft")
abline(a=0,b=1, lty=3)
c(rmse(tcdfs, ecdfs), rsquare(tcdfs, ecdfs))
}
dn <- rownames(dlf50$gof)[-(16:17)]
op <- par(mfrow=c(4,4), mar=rep(1,4), xaxt="n", yaxt="n")
for(i in dn) compranks(i)
par(op)
# so revgum, nor and rice systematically deviate from ECDF.
# RMSE is indeed better to sort by than R2


\dontrun{
dev.new(record=TRUE)
for(i in 1:10/10) distLplot(distLgof(dlf100, gofProp=i, plot=FALSE),
                            main=paste("upper", i*100, "% used for R2"))
# depending on which proportion of the data the GOF is calculated with, different
# distributions are selected to be the 5 best.
dev.off()
}
}
\keyword{univar}
\keyword{hplot}
\keyword{distribution}
