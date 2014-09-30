\name{rainbow2}
\alias{rainbow2}
\title{rainbow from blue to red}
\description{Reversed \code{\link{rainbow}} with different defaults, resulting in a color vector from blue (good) to red (bad)}
\usage{rainbow2(n=10, s=1, v=1, start=0, end=0.7, alpha=1)}
\arguments{
  \item{n}{number of colors. DEFAULT: 10}
  \item{s,v}{saturation and value as in \code{\link{rainbow}}. DEFAULT: 1}
  \item{start}{start color. DEFAULT: 0}
  \item{end}{end color. DEFAULT: 0.7}
  \item{alpha}{transparency. DEFAULT: 1)}
}

\value{A character vector of color names.}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Sept 2014}
\seealso{\code{\link{rainbow}} }
\examples{
plot(1:10, pch=16, cex=2, col=rainbow2(10))
}
\keyword{color}
\keyword{dplot}
