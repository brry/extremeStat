\name{q_evir}
\alias{q_evir}
\alias{q_evir2}
\title{GPD quantile of sample}
\description{Compute quantile of General Pareto Distribution fitted to sample
by peak over treshold (POT) method declaring treshold from truncation proportion
using the established evir (extreme values in r) package}
\usage{q_evir(x, probs, truncate, pngdev=TRUE, ...)}
\arguments{
  \item{x}{Vector with values}
  \item{probs}{Probabilities of truncated (Peak over treshold) quantile}
  \item{truncate}{Truncation percentage (proportion of sample discarded)}
  \item{pngdev}{sink \code{evir::quant} graph output to file (is removed later) instead of openening \code{\link{dev.new}}, which also is closed later. Using TRUE avoids the graphics device showing, but \code{getwd()} must be writable. DEFAULT: TRUE}
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
annMax <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6,
81.9, 60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1,
40.4, 48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)
probs2 <- c(0.8, 0.9, 0.95, 0.99, 0.999)

# Theoretically, the tails of distributions converge to GPD (General Pareto)
# Comparing the evir (extreme values in r) package with lmomco:

q_evir (x=annMax, probs=probs2, truncate=0.5)
q_evir2(x=annMax, probs=probs2, truncate=0.5)
###q_evir3(x=annMax, probs=probs2, truncate=0.5)
distLquantile(annMax, sel="gpa", truncate=0.5, probs=probs2, emp=FALSE, quiet=TRUE)
# higher estimates in this case
max(annMax)

# with bigger (random) datasets:
dlf <- distLfit(annMax)
set.seed(3) # for many seeds, quant:gpd:optim complains
rnum <- rlmomco(n=1e3, para=dlf$parameter$gpa)
#
q_evir (x=rnum, probs=probs2, truncate=0.5)
q_evir2(x=rnum, probs=probs2, truncate=0.5)
distLquantile(rnum, sel="gpa", truncate=0.5, probs=probs2, emp=FALSE, quiet=TRUE)
#
q_evir (x=rnum, probs=probs2, truncate=0.8)
q_evir2(x=rnum, probs=probs2, truncate=0.8)
distLquantile(rnum, sel="gpa", truncate=0.8, probs=probs2, emp=FALSE, quiet=TRUE)
#
# True values:
quagpa(f=probs2, para=dlf$parameter$gpa)
mt <- seq(0, 0.99, length=30)
mtq <- sapply(mt, function(t) distLquantile(rnum, sel="gpa", truncate=t,
                                         probs=probs2, emp=FALSE, quiet=TRUE))
plot(mt, mtq[5,], type="l", ylim=range(mtq, finite=TRUE), las=1)
abline(h=quagpa(f=probs2, para=dlf$parameter$gpa), col=2)
for(i in 1:5) lines(mt, mtq[i,])
# in this particular sample, distLquantile overestimates


set.seed(4)
rnum <- rlmomco(n=1e3, para=dlf$parameter$gpa)
#
q_evir (x=rnum, probs=probs2, truncate=0.5)
q_evir2(x=rnum, probs=probs2, truncate=0.5)
distLquantile(rnum, sel="gpa", truncate=0.5, probs=probs2, emp=FALSE, quiet=TRUE)
#
q_evir (x=rnum, probs=probs2, truncate=0.8)
q_evir2(x=rnum, probs=probs2, truncate=0.8)
distLquantile(rnum, sel="gpa", truncate=0.8, probs=probs2, emp=FALSE, quiet=TRUE)
#
# True values:
quagpa(f=probs2, para=dlf$parameter$gpa)
mt <- seq(0, 0.99, length=30)
mtq <- sapply(mt, function(t) distLquantile(rnum, sel="gpa", truncate=t,
                                         probs=probs2, emp=FALSE, quiet=TRUE))
plot(mt, mtq[5,], type="l", ylim=range(mtq, finite=TRUE), las=1)
abline(h=quagpa(f=probs2, para=dlf$parameter$gpa), col=2)
for(i in 1:5) lines(mt, mtq[i,])
# in this particular sample, distLquantile underestimates


data(danish)
danish <- as.numeric(danish)
#
q_evir (x=danish, probs=probs2, truncate=0.5)
q_evir2(x=danish, probs=probs2, truncate=0.5)
distLquantile(danish, sel="gpa", truncate=0.5, probs=probs2, emp=FALSE, quiet=TRUE)
#
q_evir (x=danish, probs=probs2, truncate=0.8)
q_evir2(x=danish, probs=probs2, truncate=0.8)
distLquantile(danish, sel="gpa", truncate=0.8, probs=probs2, emp=FALSE, quiet=TRUE)
#
mt <- seq(0, 0.99, length=30)
mtq <- sapply(mt, function(t) distLquantile(danish, sel="gpa", truncate=t,
                                         probs=probs2, emp=FALSE, quiet=TRUE))
plot(mt, mtq[5,], type="l", ylim=lim0(150), las=1)
abline(h=quantileMean(danish, probs=probs2), col=2)
for(i in 1:5) lines(mt, mtq[i,])
# in this particular sample, distLquantile underestimates relative to empirical quantiles

distLquantile(log10(danish), sel="gpa", plot=TRUE, probs=probs2, emp=FALSE, quiet=TRUE)
quantileMean(log10(danish), probs=probs2) # there's one extreme extreme
}
\keyword{distribution}
\keyword{robust}
\keyword{univar}
