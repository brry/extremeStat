% manually change usage section cubic meters to \U00B3
\name{plotLextreme}
\alias{plotLextreme}
\title{Plot extreme value statistics}
\usage{
plotLextreme(dlf, selection = NULL, order = FALSE, add = FALSE,
  nbest = 5, log = FALSE, xlim = NULL, ylim = NULL, las = 1,
  main = dlf$datname, xlab = "Return Period RP  [a]",
  ylab = "Discharge HQ  [m\U00B3/s]", PPcol = "black", PPpch = c(16, 3),
  PPcex = 1, distcols = berryFunctions::rainbow2(nbest), lty = 1,
  lwd = 1, pch = NA, cex = 1, n_pch = 15, legend = TRUE,
  rmse = 4, legargs = NULL, quiet = FALSE, logargs = NULL, ...)
}
\arguments{
\item{dlf}{List as returned by \code{\link{distLextreme}} or \code{\link{distLexBoot}}}

\item{selection}{Selection of distributions. Character vector with type as in
\code{\link[lmomco]{lmom2par}}. DEFAULT: NULL}

\item{order}{If selection is given, should legend and colors be ordered
by gof anyways? DEFAULT: FALSE}

\item{add}{If TRUE, plot is not called before adding lines. This lets you
add lines to an existing plot. DEFAULT: FALSE}

\item{nbest}{Number of distributions plotted, in order of goodness of fit.
Overwritten internally if selection is given. DEFAULT: 5}

\item{log}{Logarithmic x-axis? DEFAULT: FALSE}

\item{xlim}{X-axis limits. DEFAULT: xlim of plotting positions}

\item{ylim}{Y-lim. DEFAULT: from min to extended max}

\item{las}{LabelAxisStyle to orient labels, see \code{\link{par}}. DEFAULT: 1}

\item{main}{Title of plot. DEFAULT: dlf$datname}

\item{xlab}{X axis label. DEFAULT: "Return Period RP  [a]"}

\item{ylab}{Y axis label. Please note that the ubuntu pdf viewer might be
unable to display unicode superscript. DEFAULT: "Discharge HQ  [m3/s]"}

\item{PPcol}{Plotting Position point colors, vector of length two for
Weibull and Gringorton, recycled. PP are not used for
fitting distributions, but for plotting only. DEFAULT: "black"}

\item{PPpch}{point characters for plotting positions after Weibull and
Gringorton, respectively. NA to suppress in plot and legend. DEFAULT: c(16,3)}

\item{PPcex}{Character EXpansion of plotting points. DEFAULT: 1}

\item{distcols}{Color for each distribution added with \code{\link{lines}}.
Recycled, if necessary. DEFAULT: \code{\link[berryFunctions]{rainbow2}}}

\item{lty}{Line TYpe for plotted distributions. Is recycled to from a
vector of length nbest, i.e. a value for each dist. DEFAULT: 1}

\item{lwd}{Line WiDth of distribution lines. Recycled vector of length nbest.
DEFAULT: 1}

\item{pch}{Point CHaracter of points added at regular intervals.
This makes lines more distinguishable from each other.
NA to suppress. Recycled vector of length nbest. DEFAULT: NA}

\item{cex}{if pch != NA, size of points. Recycled vector of length nbest.
DEFAULT: 1}

\item{n_pch}{Number of points spread evenly along the line.
Recycled vector of length nbest. DEFAULT: 15}

\item{legend}{Logical. Add a legend? DEFAULT: TRUE}

\item{rmse}{Integer. If rmse > 0, RMSE values are added to legend.
They are rounded to \code{rmse} digits. DEFAULT: 4}

\item{legargs}{list of arguments passed to \code{\link{legend}} except for
legend, col, pch, lwd, lty. DEFAULT: NULL}

\item{quiet}{Suppress notes? DEFAULT: FALSE}

\item{logargs}{list of arguments passed to \code{berryFunctions::\link{logAxis}}.}

\item{\dots}{Further arguments passed to \code{\link{plot}} like yaxt="n", ...}
}
\value{
invisible dlf object, see \code{\link{printL}}
}
\description{
Plots distributions fitted by L-moments and adds plotting positions by Weibull and Gringorton.
This is an auxiliary graphing function to \code{\link{distLextreme}}
}
\examples{
#see
?distLextreme

}
\seealso{
\code{\link{distLextreme}}, \code{\link{plotLfit}}
}
\author{
Berry Boessenkool, \email{berry-b@gmx.de}, March 2015, updated heavily Aug 2015
}
\keyword{distribution}
\keyword{dplot}
\keyword{hplot}
