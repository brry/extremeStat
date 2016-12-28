

#' annual discharge maxima (streamflow)
#' 
#' Annual discharge maxima of a stream in Austria called Griesler or Fuschler
#' Ache, at the measurement station (gauge) near St. Lorenz, catchment area ca
#' 100 km^2. Extracted from the time series 1976-2010 with a resolution of 15
#' Minutes.
#' 
#' 
#' @name annMax
#' @docType data
#' @format num [1:35] 61.5 77 37 69.3 75.6 74.9 43.7 50.8 55.6 84.1 ...
#' @source Hydrographische Dienste Oberoesterreich und Salzburg, analysed by
#' package author (\email{berry-b@gmx.de})
#' @keywords datasets
#' @examples
#' 
#' data(annMax)
#' str(annMax)
#' str(annMax)
#' plot(1976:2010, annMax, type="l", las=1, main="annMax dataset from Austria")
#' # Moving Average with different window widths:
#' berryFunctions::movAvLines(annMax, x=1976:2010, lwd=3, alpha=0.7)
#' 
NULL




#' Extreme value statistics on a linear scale
#' 
#' Fit (via L moments), plot (on a linear scale) and compare (by goodness of fit)
#' several (extreme value) distributions.
#' Compute high quantiles even in small samples and estimate extrema at given return periods.\cr
#' Open the \href{https://cran.r-project.org/package=extremeStat/vignettes/extremeStat.html}{Vignette} 
#' for an introduction to the package: \code{vignette("extremeStat")}\cr
#' This package heavily relies on and thankfully acknowledges the package \code{lmomco} by WH Asquith.
#' 
#' @details
#' 
#' The main functions in the extremeStat package are:
#' \tabular{lll}{
#' \code{\link{distLweights}} -> \tab \code{\link{distLgof}} \tab -> \code{\link{plotLgof}} \cr
#'                               \tab \code{\link{distLfit}} \tab -> \code{\link{plotLfit}} \cr
#' \code{\link{q_gpd}} + \code{\link{q_weighted}} -> \tab \code{\link{distLquantile}} \tab -> \code{\link{plotLquantile}} \cr
#'                               \tab \code{\link{distLextreme}} \tab -> \code{\link{plotLextreme}} \cr
#'                               \tab \code{\link{distLexBoot}} \tab \cr
#' }
#' 
#' The common object to share between functions is a list that can be printed with 
#' \code{\link{printL}(dlf)} and contains the following elements:
#' \tabular{ll}{
#' \code{dat}          \tab numeric vector with (extreme) values, 
#'                          with all NAs and values below threshold removed\cr
#' \code{dat_full}     \tab original input data complete with NAs\cr
#' \code{datname}      \tab character string for main, xlab etc \cr
#' \code{parameter}    \tab list (usually of length 17 if \code{speed=TRUE} in 
#'                          \code{\link{distLfit}}) 
#'                          with parameters of each distribution\cr
#' \code{gof}          \tab dataframe with 'Goodness of Fit' measures, sorted by 
#'                          RMSE of theoretical and empirical cumulated density,
#'                          created by \code{\link{distLgof}}\cr
#' \code{distnames}    \tab character vector with selected distribution names\cr
#' \code{distcols}     \tab colors for distnames (for plotting). If not given manually,
#'                          determined by \code{berryFunctions::\link{rainbow2}}\cr
#' \code{distselector} \tab character string with function name creating
#'                          the selection\cr
#' \code{truncate, threshold} \tab Truncation percentage and threshold value, 
#'                          relevant for \code{\link{distLquantile}}\cr
#'                       
#' optionally:         \tab \cr
#' \code{returnlev, }  \tab dataframe with values of distributions for given 
#'                          return periods (\code{RPs}). This element is only added
#'                          in \code{\link{distLextreme}}\cr
#' \code{RPweibull, RPgringorton} \tab Return periods according to plotting positions, 
#'                          added in \code{\link{plotLextreme}}\cr
#' \code{quant}        \tab Quantile estimates from \code{\link{distLquantile}}\cr
#' \code{exBootRPs, qexBootSim, exBootCI} \tab objects from \code{\link{distLexBoot}}\cr

#' }
#' 
#' Plotting positions are not used for fitting distributions, but for plotting only.
#' The ranks of ascendingly sorted extreme values are used to 
#' compute the probability of non-exceedence Pn:\cr
#' \code{Pn_w <-  Rank      /(n+1)       # Weibull}\cr
#' \code{Pn_g <- (Rank-0.44)/(n+0.12)    # Gringorton (taken from lmom:::evplot.default)}\cr
#' Finally: RP = Returnperiod = recurrence interval = 1/P_exceedence = 1/(1-P_nonexc.), thus:\cr
#' \code{RPweibull = 1/(1-Pn_w)} and analogous for gringorton.\cr
#'

#' This package was motivated by my need to compare the fits of several distributions to data.
#' It was originally triggered by a flood estimation assignment we had in class 2012,
#' and it bothered me that we just assumed the gumbel distribution would fit the data fine.
#' 
#' @name extremeStat
#' @aliases extremeStat-package extremeStat
#' @docType package
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014-2016

#' @seealso
#' If you are looking for more detailed (uncertainty) analysis, eg confidence intervals,
#' check out the package \code{extRemes}, especially the function \code{\link[extRemes]{fevd}}.
#' \url{https://cran.r-project.org/package=extRemes}\cr
#' Intro slides: \url{http://sites.lsa.umich.edu/eva2015/wp-content/uploads/sites/44/2015/06/Intro2EVT.pdf}\cr
#' Parameter fitting and distribution functions: \url{https://cran.r-project.org/package=lmomco}\cr
#' Distributions: \url{https://www.rmetrics.org/files/Meielisalp2009/Presentations/Scott.pdf}
#' and: \url{https://cran.r-project.org/view=Distributions} \cr
#' R in Hydrology: \url{http://abouthydrology.blogspot.de/2012/08/r-resources-for-hydrologists.html}\cr

#' @keywords package documentation
#' @importFrom grDevices extendrange
#' @importFrom graphics abline axis box hist legend lines par plot points text
#' @importFrom methods slotNames
#' @importFrom stats approx ecdf ks.test median quantile
#' @importFrom utils head tail

#' @examples
#' 
#' data(annMax) # annual discharge maxima from a stream in Austria
#' plot(annMax, type="l")
#' dle <- distLextreme(annMax)
#' dle$returnlev
#' 
NULL



