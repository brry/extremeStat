\name{rmse}
\alias{rmse}
\title{Root Mean Squared Error}
\description{Root Mean Squared Error}
\usage{rmse(a,b)}
\arguments{
  \item{a}{Vector with values.}
  \item{b}{Another vector of the same length.}
}
\details{Formula used: \code{sqrt( sum((a-b)^2)/length(b) )}}
\value{Numeric.}
\references{\url{http://en.wikipedia.org/wiki/Mean_squared_error}}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, 2014}
\seealso{\code{\link{rsquare}} }
\examples{
x <- rnorm(20)
y <- 2*x + rnorm(20)
plot(x,y)
yp <- predict(lm(y~x))
plot(y, yp)
abline(a=0,b=1)
rmse(y,yp)
}
\keyword{univar}
