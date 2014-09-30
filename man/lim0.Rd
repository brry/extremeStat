\name{lim0}
\alias{lim0}
\title{axis limits with one end at zero}
\description{
Calculates the range needed for ylim or xlim in plot, so that axis starts at zero and is extended by 4\% at the other end}
\usage{lim0(x, f=1/27, curtail=TRUE)}
\arguments{
  \item{x}{Numeric. Vector with values}
  \item{f}{Numeric. Extension factor. DEFAULT: 0.04 as in extendrange used eg. by \code{\link{curve}}}
  \item{curtail}{Logical. Should the range returned be trimmed by 4\%? That way, plotting doesn't need the default \code{\link{par}} xaxs or yaxs changed. DEFAULT: TRUE}
}
\value{Vector with two values}
\references{\code{methods(plot)}, \code{\link[graphics]{plot.default}}.
actually, I found extendrange via \code{plot.function} in \code{curve}
}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, 6.6.2013 }
\note{Needs more testing, please get back to me if it doesn't do what it's supposed to do!}
\seealso{The \code{\link{extendrange}()} utility in package \pkg{grDevices}}
\examples{
# basic idea:
val <- c(3.2, 1.8, 4.5, 2.8, 0.1, 2.9) # just some numbers
plot(val, ylim=lim0(val) ) # you don't even have to set yaxs="i" ;-)

# "normal" plot:
plot(val)
par("usr")  # -0.076  4.676

# if y-axis is not allowed to go below 0, and we're too lazy to set yaxs="i":
plot(val, ylim=lim0(val) )
round( par("usr")  , digits=5) # 0.00000 4.66296

# with 0.04 extension as claimed by help page (1/27 in source code = 0.037):
plot(val, ylim=lim0(val, f=0.04) )
round( par("usr")  , digits=5) # zero is not included on axis anymore

b <- -val
plot(b)
plot(b, ylim=lim0(b) ) # works with only negative values as well
}
\keyword{ dplot }
