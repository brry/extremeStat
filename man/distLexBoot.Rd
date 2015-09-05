\name{distLexBoot}
\alias{distLexBoot}
\title{Uncertainty intervals for return periods}
\description{calculates and plots bootstrap uncertainty intervals for \code{\link{distLextremePlot}}.}
\usage{distLexBoot(dlf, nbest=3, selection=NULL, truncate=0, 
n=100, prop=0.8, returnall=FALSE, conf.lev=0.95, RPs=NULL, 
plot=TRUE, add=FALSE, log=TRUE, progbars=TRUE, ... )}
\arguments{
  \item{dlf}{dlf object, as returned by \code{\link{distLextreme}}, is passed to \code{\link{distLextremePlot}}.}
  \item{nbest}{Number of best fitted distribution functions in dlf for which bootstrapping is to be done. Overriden by \code{selection}. DEFAULT: 3}
  \item{selection}{Character vector with distribution function names to be used. Suggested to keep this low. DEFAULT: NULL}
  \item{truncate}{Truncation of subsamples, see \code{\link{distLquantile}}. DEFAULT: 0}
  \item{n}{Number of subsamples to be processed (computing time increases extraordinarily). DEFAULT: 100}
  \item{prop}{Proportion of sample to be used in each run. DEFAULT: 0.8}
  \item{returnall}{Return all simulations, instead of the aggregate confidence level? DEFAULT: FALSE}
  \item{conf.lev}{Confidence level (Proportion of subsamples within 'confidence interval'). Quantiles extracted from this value are passed to \code{\link[berryFunctions]{quantileMean}}. DEFAULT: 0.95}
  \item{RPs}{Return Period vector, by default calculated internally based on log. DEFAULT: NULL}
  \item{plot}{Plot results via \code{\link{distLextremePlot}}? DEFAULT: TRUE}
  \item{add}{Add to existing plot? DEFAULT: FALSE}
  \item{log}{Plot on a logarithmic axis. DEFAULT: TRUE}
  \item{progbars}{Show progress bar for Monte Carlo simulation? DEFAULT: TRUE}
  \item{\dots}{Further arguments passed to \code{\link{distLextremePlot}}}
}
\details{Has not been thoroughly tested yet. 
Bootstrapping defaults can probabply be improved.}
\value{a list with for each selection, a matrix with confidence intervals for RPs, or if returnall=TRUE, all the simulation results}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2015}
\seealso{\code{\link{distLextreme}} }
\examples{
data(annMax)
dlf <- distLextreme(annMax, log=TRUE, selection=c("wak","gum","gev","nor"))
dleB <- distLexBoot(dlf, nbest=4, conf.lev=0.5)
}
\keyword{bootstrap}
\keyword{montecarlo}
\keyword{hplot}
\keyword{dplot}
\keyword{distribution}
\keyword{ts}