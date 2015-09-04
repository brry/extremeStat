\name{annMax}
\alias{annMax}
\docType{data}
\title{annual discharge maxima (streamflow)}
\description{Annual discharge maxima of a stream in Austria called Griesler or Fuschler Ache, at the measurement station (gauge) near St. Lorenz, catchment area ca 100 km^2.
Extracted from the time series 1976-2010 with a resolution of 15 Minutes.}
\usage{data("annMax")}
\format{num [1:35] 61.5 77 37 69.3 75.6 74.9 43.7 50.8 55.6 84.1 ...}
\source{Hydrographische Dienste Oberoesterreich und Salzburg, analysed by package author (\email{berry-b@gmx.de})}
\examples{
data(annMax)
str(annMax)
str(annMax)
plot(1976:2010, annMax, type="l", las=1, main="annMax dataset from Austria")
# Moving Average with different window widths:
berryFunctions::movAvLines(annMax, x=1976:2010, lwd=3, alpha=0.7)
}
\keyword{datasets}
