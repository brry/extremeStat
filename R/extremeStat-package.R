

#' Extreme value statistics on a linear scale
#' 
#' Fit (via L moments), plot (on a linear scale) and compare (by goodness of fit)
#' several (extreme value) distributions.
#' Compute high quantiles even in small samples and estimate extrema at given return periods.\cr
#' Open the \href{https://cran.r-project.org/package=extremeStat/vignettes/extremeStat.html}{Vignette} 
#' for an introduction to the package: \code{vignette("extremeStat")}\cr
#' This package heavily relies on and thankfully acknowledges the package \code{lmomco} by WH Asquith.
#' 
#' @section Package overview:
#' 
#' The main functions in the extremeStat package are:
#' \tabular{ll}{
#' \code{\link{distLweights}} \tab -> \code{\link{plotLweights}}    \cr
#' \code{\link{distLfit}}     \tab -> \code{\link{plotLfit}}     \cr
#' \code{\link{q_gpd}} + \code{\link{q_weighted}} -> \code{\link{distLquantile}} \tab -> \code{\link{plotLquantile}} \cr
#' \code{\link{distLextreme}} \tab -> \code{\link{plotLextreme}} \cr
#' \code{\link{distLexBoot}}  \tab \cr
#' }
#' They create and modify a list object printed by (and documented in) \code{\link{printL}}.
#' 
#' @name extremeStat
#' @aliases extremeStat-package extremeStat
#' @docType package
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014-2016
#' 
#' @seealso
#' If you are looking for more detailed (uncertainty) analysis, eg confidence intervals,
#' check out the package \code{extRemes}, especially the function \code{\link[extRemes]{fevd}}.
#' \url{https://cran.r-project.org/package=extRemes}\cr
#' Intro slides: \url{http://sites.lsa.umich.edu/eva2015/wp-content/uploads/sites/44/2015/06/Intro2EVT.pdf}\cr
#' Parameter fitting and distribution functions: \url{https://cran.r-project.org/package=lmomco}\cr
#' Distributions: \url{https://www.rmetrics.org/files/Meielisalp2009/Presentations/Scott.pdf}
#' and: \url{https://cran.r-project.org/view=Distributions} \cr
#' R in Hydrology: \url{http://abouthydrology.blogspot.de/2012/08/r-resources-for-hydrologists.html}\cr
#' 
#' @keywords package documentation
#' @importFrom grDevices extendrange
#' @importFrom graphics abline axis box hist legend lines par plot points text
#' @importFrom methods slotNames
#' @importFrom stats approx ecdf ks.test median quantile
#' @importFrom utils head tail
#' 
#' @examples
#' data(annMax) # annual discharge maxima from a stream in Austria
#' plot(annMax, type="l")
#' dle <- distLextreme(annMax)
#' dle$returnlev
#' 
NULL




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



#' distribution weights
#' 
#' Weights for weighted average as in the submission of revisions for the paper
#' \url{http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183/}
#' 
#' 
#' @name weightp
#' @docType data
#' @format named num [1:17] 
#' @source See paper revisions (not yet online at moment of extremeStat update) (\email{berry-b@gmx.de})
#' @keywords datasets
#' @examples
#' 
#' data(weightp)
#' data.frame(weightp)
#' barplot(weightp, horiz=TRUE, las=1)
#' stopifnot(   all.equal(sum(weightp), 1)   )
#' 
if(FALSE){
weightp1 <- read.table(text="
wei    0.12915523
pe3    0.12645899
gpa    0.12070367
wak    0.11883859
gno    0.09569464
gum    0.08812947
exp    0.07123161
gev    0.06609407
lap    0.06180443
ray    0.06040467
gam    0.04881928
glo    0.01266536
rice   0.00000000
nor    0.00000000
ln3    0.00000000
revgum 0.00000000
kap    0.00000000")
weightp <- weightp1[,2]
names(weightp) <- weightp1[,1]
rm(weightp1)
devtools::use_data(weightp)
}
"weightp"

