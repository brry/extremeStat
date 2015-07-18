\name{distLfit}
\alias{distLfit}
\title{Fit distributions via linear moments}
\description{Fit several distributions via linear moments, plot histogram and distribution densities \emph{or} ecdf with cumulated probability.}
\usage{distLfit(dat, datname, speed=TRUE, ks=FALSE, selection=NULL,
gofProp=1, truncate=0, gofComp=FALSE, progbars=length(dat)>200, time=TRUE,
plot=TRUE, cdf=FALSE, legargs=NULL, histargs=NULL, quiet=FALSE, ... )}
\arguments{
  \item{dat}{Vector with values}
  \item{datname}{Character string for main, xlab etc. DEFAULT: internal \code{substitute(dat)}}
  \item{speed}{If TRUE, several distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}. DEFAULT: TRUE}
  \item{ks}{Include ks.test results in \code{dlf$gof}? Computing is much faster when FALSE. DEFAULT: FALSE}
  \item{selection}{Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL}
  \item{gofProp}{Upper proportion (0:1) of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail. DEFAULT: 1}
  \item{truncate}{Number between 0 and 1. Censored \code{\link{distLquantile}}: fit to highest values only (truncate lower proportion of x). Probabilities are adjusted accordingly. DEFAULT: 0}
  \item{gofComp}{If TRUE, plots a comparison of the ranks of different GOF-methods and sets plot to FALSE. DEFAULT: FALSE}
  \item{progbars}{Show progress bars for each loop? DEFAULT: TRUE if n > 200}
  \item{time}{\code{\link{message}} execution time? DEFAULT: TRUE}
  \item{plot}{Should a histogram with densities be plotted? DEFAULT: TRUE}
  \item{cdf}{If TRUE, plot cumulated DF instead of probability density. DEFAULT: FALSE}
  \item{legargs}{List of arguments passed to \code{\link{legend}} except for legend and col. DEFAULT: NULL}
  \item{histargs}{List of arguments passed to \code{\link{hist}} except for x, breaks, col, xlim, freq. DEFAULT: NULL}
  \item{quiet}{Suppress notes? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link{distLplot}} if they are accepted there, else passed to \code{\link{lines}}, like lty, type, pch, ...}
}
\details{fits parameters via \code{\link[lmomco]{lmom2par}} in the package \code{lmomco}}
\value{List as explained in \code{\link{extremeStat}}.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{distLgof}}, \code{\link{distLplot}}. \code{\link[extRemes]{fevd}} in the package \code{extRemes},
\code{\link[MASS]{fitdistr}} in the package \code{MASS}.}
\examples{
testdat <- c(61.5, 77.0, 37.0, 69.3, 75.6, 74.9, 43.7, 50.8, 55.6, 84.1, 43.6, 81.9,
60.1, 72.4, 61.6, 94.8, 82.6, 57.2,  63.1, 73.8, 51.3, 93.6, 56.9, 52.1, 40.4,
48.9, 113.6, 35.4, 40.1, 89.6, 47.8, 57.6, 38.9, 69.7, 110.8)

# basic usage on real data (annual discharge maxima in Austria)
dlf <- distLfit(testdat)
str(dlf, max.lev=2)
distLprint(dlf)

# arguments that can be passed:
distLfit(testdat, lty=2, col=3, legargs=list(lwd=3), main="booh!")
dlf_b <- distLfit(rbeta(1e3, 5, 2), nbest=10, legargs=c(x="left"))
distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"))
distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"), order=TRUE)

# Goodness of Fit is computed by RMSE, see first example of ?distLgof

# logarithmic axes:
set.seed(1)
y <- 10^rnorm(1e3, mean=2, sd=0.3) # if you use 1e4, distLgof will be much slower
hist(y, breaks=20)
hist(log10(y))
dlf <- distLfit(log10(y), breaks=50)
distLplot(dlf, breaks=50, log=TRUE)

\dontrun{
# this takes a while, as it tries to fit all 30 distributions:
d_all <- distLfit(testdat, gofProp=1, speed=FALSE, plot=FALSE) # 35 sec
distLprint(d_all)
distLplot(d_all, nbest=22, coldist=grey(1:22/29), xlim=c(20,140))
distLplot(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
d_all$gof
}
}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
