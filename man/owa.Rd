\name{owa}
\alias{owa}
\title{overwrite argument lists passed to functions}
\description{Can be used in functions that give defaults for several argument lists that are each passed to different functions.
Some of the defaults can be overwritten, some should be left unchanged, some are
aditionally specified by users. owa combines everything accordingly.}
\usage{owa(d, a, u=NA)}
\arguments{
  \item{d}{Default arguments}
  \item{a}{Arguments specified by user}
  \item{u}{names of Unchangeable arguments (that can not be overwritten). DEFAULT: NA}
}
\details{None}
\value{Always a list, disregarding list/vector mode of input}
\author{Berry Boessenkool, \email{berry-b@gmx.de}, Early 2014}
\references{
   \url{http://stackoverflow.com/questions/3057341}\cr
   \url{http://stackoverflow.com/questions/5890576}\cr
   \url{http://stackoverflow.com/questions/4124900}\cr
   \url{http://stackoverflow.com/questions/16774946}\cr}
\examples{
d <- list(col="gray", lty=1, wachs=1:6, bb=42)
a <- list(lty=3, lwd=2, wachs="A", bb=16)

owa(d, a, u=c("bb", "wachs"))
owa(d, NULL, u=c("bb", "wachs")) # NULL is a good default for argument lists
owa(d, c(HH=2, BBB=3) ) # vectors and lists all converted to lists


# An example of how to code (as in funnelPlot or mReg):
# Why we want to do this:
testfun <- function(...) {plot(7:9, ...) ; legend("top", "Text hier", ...)}
testfun()
\dontrun{
testfun(type="o") # legend doesn't have the argument 'type'!
}

# How we do this:
testfun <- function(data=7:9, legarg=NULL, plotarg=NULL)
   {
   # defaults for plot and legend:
   plot_def <- list(x=0.5*data, col="red", cex=2, lty=2, type="o")
   leg_def <- list(x="top", lty=2, legend="Default text here")
   # combine defaults and user specified into final argument list
   plot_fin <- owa(d=plot_def, a=plotarg, u=c("col", "lty"))
   leg_fin <- owa(d=leg_def, a=legarg, u="lty")
   # Execute single functions that each have their own arguments:
   do.call(  plot, args=plot_fin)
   do.call(legend, args=leg_fin)
   }

testfun()
testfun(plotarg=list(type="l", col="blue") )
# color is silently ignored, as it is defined as unchangeable
testfun(plotarg=list(type="l"), legarg=list(col="blue", pch=16) )

}
\keyword{programming}
