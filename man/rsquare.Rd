\name{rsquare}
\alias{rsquare}
\title{R squared}
\description{R squared (coefficient of determination)}
\usage{rsquare(a,b, quiet=FALSE)}
\arguments{
  \item{a}{Vector with values.}
  \item{b}{Another vector of the same length.}
  \item{quiet}{Should NA-removal warnings be suppressed? Helpful within functions. DEFAULT: FALSE}
}
\details{Formula used: \code{\link{cor}(a,b)^2}}
\value{Numeric.}
\references{\url{http://en.wikipedia.org/wiki/R-squared}}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, 2014}
\note{Using cor is much faster than using\cr
\code{  aa <-  a-mean(a);
  bb <-  b-mean(b);
  sum(aa*bb)^2/sum(aa^2)/sum(bb^2)}}
\seealso{\code{\link{rmse}}, \code{\link{cor}}, \code{\link{lm}}}
\examples{
x <- rnorm(20)
y <- 2*x + rnorm(20)
plot(x,y)
rsquare(x,y)

r2 <- sapply(1:10000, function(i){
   x <- rnorm(20);  y <- 2*x + rnorm(20);  rsquare(x,y) })
hist(r2, breaks=70, col=5,
main= "10'000 times   x <- rnorm(20);  y <- 2*x + rnorm(20);  rsquare(x,y)")

}
\keyword{univar}
