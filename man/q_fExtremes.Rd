\name{q_fExtremes}
\alias{q_fExtremes}
\title{GPD quantile of sample}
\description{Compute quantile of General Pareto Distribution fitted to sample
by peak over treshold (POT) method declaring treshold from truncation proportion
using the fExtremes package (Extreme Financial Market Data)}
\usage{q_fExtremes(x, probs, truncate, undertruncNA=TRUE, method="pwm", quiet=FALSE, ...)}
\arguments{
  \item{x}{Vector with values. NAs are silently ignored}
  \item{probs}{Probabilities of truncated (Peak over treshold) quantile}
  \item{truncate}{Truncation percentage (proportion of sample discarded)}
  \item{undertruncNA}{Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE}
  \item{method}{Method in \code{\link[fExtremes]{gpdFit}}, "mle" (maximum likelihood estimation) or "pwm" (probability-weighted moments). Note this is called type in the original function! DEFAULT: "pwm"}
  \item{quiet}{Should messages be suppressed? DEFAULT: FALSE}
  \item{\dots}{Further arguments passed to \code{\link[fExtremes]{gpdFit}}, except for \code{type}, as it is already specified by \code{method}.}
}
\value{Vector of quantile estimates for each value of probs}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Feb 2016}
\seealso{\code{\link{distLquantile}}, \code{\link{q_evir}}, \code{\link[fExtremes]{gpdFit}}}
\examples{
library(fExtremes)
data(annMax)
p <- c(0.5, 0.8, 0.9, 0.99, 0.999)
d <- q_fExtremes(annMax, truncate=0,   p=p); d
d <- q_fExtremes(annMax, truncate=0.6, p=p); d
d <- q_fExtremes(annMax, truncate=0.6, p=p, under=FALSE); d
d <- q_fExtremes(annMax, truncate=1,   p=p); d
d <- q_fExtremes(annMax, truncate=1,   p=p, under=FALSE); d
}
}
\keyword{distribution}
\keyword{robust}
\keyword{univar}
