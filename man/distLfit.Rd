\name{distLfit}
\alias{distLfit}
\title{Fit distributions via linear moments}
\description{Fit several distributions via linear moments, plot histogram and distribution densities \emph{or} ecdf with cumulated probability.}
\usage{distLfit(dat, datname, speed=TRUE, ks=FALSE, selection=NULL,
gofProp=1, truncate=0, gofComp=FALSE, progbars=length(dat)>200, time=TRUE,
plot=TRUE, cdf=FALSE, legargs=NULL, histargs=NULL, quiet=FALSE, quietss=quiet, ... )}
\arguments{
  \item{dat}{Vector with values}
  \item{datname}{Character string for main, xlab etc. DEFAULT: internal \code{substitute(dat)}}
  \item{speed}{If TRUE, several distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list}()}. DEFAULT: TRUE}
  \item{ks}{Include ks.test results in \code{dlf$gof}? Computing is much faster when FALSE. DEFAULT: FALSE}
  \item{selection}{Selection of distributions. Character vector with types as in \code{\link[lmomco]{lmom2par}}. Overrides speed. DEFAULT: NULL}
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
  \item{quietss}{Suppress sample size notes? DEFAULT: quiet}
  \item{\dots}{Further arguments passed to \code{\link{distLplot}} if they are accepted there, else passed to \code{\link{lines}}, like lty, type, pch, ...}
}
\details{fits parameters via \code{\link[lmomco]{lmom2par}} in the package \code{lmomco}}
\value{List as explained in \code{\link{extremeStat}}.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{distLgof}}, \code{\link{distLplot}}. \code{\link[extRemes]{fevd}} in the package \code{extRemes},
\code{\link[MASS]{fitdistr}} in the package \code{MASS}.}
\examples{
data(annMax)
# basic usage on real data (annual discharge maxima in Austria)
dlf <- distLfit(annMax)
str(dlf, max.lev=2)
distLprint(dlf)

# arguments that can be passed:
distLfit(annMax, lty=2, col=3, legargs=list(lwd=3), main="booh!")
dlf_b <- distLfit(rbeta(1e3, 5, 2), nbest=10, legargs=c(x="left"))
distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"))
distLplot(dlf_b, selection=c("gpa", "glo", "gev", "wak"), order=TRUE)
distLplot(dlf_b, coldist=c("orange",3:6), lty=1:3) # lty is recycled
distLplot(dlf_b, cdf=TRUE)
distLplot(dlf_b, cdf=TRUE, histargs=list(do.points=FALSE), sel="nor")


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
d_all <- distLfit(annMax, gofProp=1, speed=FALSE, plot=FALSE) # 35 sec
distLprint(d_all)
distLplot(d_all, nbest=22, coldist=grey(1:22/29), xlim=c(20,140))
distLplot(d_all, nbest=22, histargs=list(ylim=c(0,0.04)), xlim=c(20,140))
d_all$gof
}
}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
