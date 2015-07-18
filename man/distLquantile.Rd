\name{distLquantile}
\alias{distLquantile}
\title{distribution quantiles}
\description{Parametric quantiles of distributions fitted to a sample.}
\usage{distLquantile(x=NULL, probs=0:4/4, truncate=0, selection=NULL, type=NULL,
dlf=NULL, order=TRUE, plot=FALSE, cdf=FALSE, lines=TRUE, linargs=NULL,
empirical=TRUE, weighted=TRUE, quiet=FALSE, trans=quiet, ... )}
\arguments{
  \item{x}{Sample for which parametrical quantiles are to be calculated. If it is NULL (the default), \code{dat} from \code{dlf} is used. DEFAULT: NULL}
  \item{probs}{Numeric vector of probabilities with values in [0,1]. DEFAULT: 0:4/4}
  \item{truncate}{Number between 0 and 1. Censored quantile: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly. DEFAULT: 0}
  \item{selection}{Distribution type, eg. "gev" or "wak", see \code{\link[lmomco]{dist.list} in lmomco}. Can be a vector. If NULL (the default), all types present in dlf$parameter are used. DEFAULT: NULL}
  \item{type}{Overrides selection. Kept for backwards compatibility. (Changed to selection for internal consistency). DEFAULT: NULL}
  \item{dlf}{dlf object described in \code{\link{extremeStat}}. Ignored if x is given, or if truncate / selection do not match. Use this to save computing time for large datasets where you already have dlf. DEFAULT: NULL}
  \item{order}{Sort results by GOF? If TRUE (the default) and length(type)>1, the output is ordered by dlf$gof, else by order of appearance in dlf$parameter. DEFAULT: TRUE}
  \item{plot}{Should \code{\link{distLplot}} be called? DEFAULT: FALSE}
  \item{cdf}{If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE}
  \item{lines}{Should vertical lines marking the quantiles be added? DEFAULT: TRUE}
  \item{linargs}{Arguments passed to \code{\link{lines}}. DEFAULT: NULL}
  \item{empirical}{Add vertical line for empirical \code{\link{quantileMean}} (and include the result in the output matrix)? DEFAULT: TRUE}
  \item{weighted}{Include weighted averages across distribution functions to the output? DEFAULT: TRUE}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{trans}{Suppress note about transposing? DEFAULT: quiet}
  \item{\dots}{Arguments passed to \code{\link{distLfit}}.}
}
\details{
Very high quantiles (99\% and higher) need large sample sizes for \code{\link{quantile}} to yield a robust estimate.
Theoretically, at least 1/(1-probs) values must be present, e.g. 10'000 for Q99.99\%.
With smaller sample sizes (eg n=35), they underestimate the actual (but unknown) quantile.
Parametric quantiles need only small sample sizes. They don't have a systematical underestimation bias, but have higher variability.
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
# Annual Discharge Maxima (streamflow)
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6,
81.9, 60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1,
40.4, 48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)
dlf <- distLfit(annMax)

distLquantile(annMax)
distLquantile(annMax, probs=0.95, plot=TRUE, linargs=c(lwd=2), nbest=5, breaks=10)
# Parametric 95% quantile estimates range from 92 to 111!
# But the best fitting distributions all lie aroud 103.

# weighted distribution quantiles are calculated by different weighting schemes:
distLgofPlot(dlf, ranks=FALSE, weights=TRUE)

# The default is to fit parameters. If speed is important in big datasets
# and parameters are already available, pass them via dlf:
distLquantile(dlf=dlf, probs=0:10/10, selection=c("wak","gev","kap"), order=FALSE)



# censored (truncated, trimmed) quantile:
qwak <- distLquantile(annMax, sel="wak", prob=0.95, plot=TRUE, ylim=c(0,0.06), emp=FALSE)
qwak2 <-distLquantile(annMax, sel="wak", prob=0.95, truncate=0.6, plot=TRUE,
                     add=TRUE, col=6, breaks=10, coldist="blue", empirical=FALSE)

# Simulation of truncation effect
library(lmomco)
#set.seed(42)
rnum <- rlmomco(n=1e3, para=dlf$parameter$gev)
myprobs <- c(0.9, 0.95, 0.99, 0.999)
mytrunc <- seq(0, 0.9, length.out=20)
trunceffect <- sapply(mytrunc, function(mt) distLquantile(rnum, selection="gev",
                             probs=myprobs, truncate=mt, plot=FALSE, trans=TRUE,
                             progbars=FALSE, empirical=FALSE, weight=FALSE))
# If more values are truncated, the function runs faster

op <- par(mfrow=c(2,1), mar=c(2,4.5,2,0.5), cex.main=1)
distLquantile(rnum, sel="gev", probs=myprobs, emp=FALSE, progbars=FALSE,
              ylab="", xlab="", plot=TRUE)
distLquantile(rnum, sel="gev", probs=myprobs, emp=FALSE, progbars=FALSE,
              truncate=0.3, add=TRUE, coldist=4, plot=TRUE)
legend("right", c("fitted GEV", "fitted with truncate=0.3"), lty=1, col=c(2,4),
       bg="white")
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


# Developmental testing:
distLquantile(annMax, plot=FALSE, selection="wak", empirical=FALSE)
distLquantile(annMax, plot=TRUE, selection="wak", empirical=FALSE, breaks=10)

distLquantile(rexp(199), sel=c("wak", "gpa"), truncate=0.8, probs=c(0.7, 0.8, 0.9))
# distLquantile(rexp(199), truncate=0.8, probs=0.7) # correctly throws an error
distLquantile(rexp(199), selection=c("wak", "gpa"))
distLquantile(rexp(199), selection="gpa")
distLquantile(rexp(4))
distLquantile(rexp(4), selection="gpa") # should not be an Error anymore
dist.list()
# distLquantile(rexp(199), selection=1:5, emp=FALSE) # index is a bad idea anyways
# distLquantile(rexp(199), selection=-3)

}
\keyword{distribution}
\keyword{robust}
\keyword{univar}
