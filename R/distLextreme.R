#' Extreme value stats
#' 
#' Extreme value statistics for flood risk estimation.
#' Input: vector with annual discharge maxima.
#' Output: discharge estimates for given return periods, 
#' parameters of several distributions (fit based on linear moments), 
#' quality of fits, plot with linear axis (dists + plotting positions by Weibull and Gringorton).
#'
#' @details \code{\link{distLextremePlot}} adds weibull and gringorton plotting positions
#' to the distribution lines, which are estimated from the linear moments of the data itself.\cr
#' I personally believe that if you have, say, 35 values in \code{dat},
#' the highest return period should be around 36 years (Weibull) and not 60 (Gringorton).\cr
#' The plotting positions don't affect the distribution parameter estimation, so this dispute is not really important.
#' But if you care, go ahead and google "weibull vs gringorton plotting positions".
#'
#' @return List as explained in \code{\link{extremeStat}}. The added element is
#'         \code{returnlev}, a data.frame with the return level (discharge) for all given RPs and for each distribution.
#' @note This function replaces \code{berryFunctions::extremeStatLmom}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012 (first draft) - 2014 & 2015 (main updates)
#' @seealso \code{\link{distLfit}}. \code{\link{distLexBoot}} for confidence
#'          interval from Bootstrapping. \code{\link[extRemes]{fevd}} in the package \code{extRemes}.
#' @references \url{http://RclickHandbuch.wordpress.com} Chapter 15 (German)\cr
#'             Christoph Mudersbach: Untersuchungen zur Ermittlung von hydrologischen
#'             Bemessungsgroessen mit Verfahren der instationaeren Extremwertstatistik
#' @keywords hplot dplot distribution ts
#' @export
#' @importFrom berryFunctions owa
#' 
#' @examples
#' 
#' # Basic examples
#' # Advanced options
#' # weighted mean based on Goodness of fit (GOF)
#' # Effect of data proportion used to estimate GOF
#' # compare extremeStat with other packages
#' 
#' library(lmomco)
#' library(berryFunctions)
#' 
#' data(annMax) # annual streamflow maxima in river in Austria
#' 
#' # Basic examples ---------------------------------------------------------------
#' dle <- distLextreme(annMax, log=TRUE)
#'
#' # Object structure:
#' str(dle, max.lev=2)
#' distLprint(dle)
#' 
#' # discharge levels for default return periods:
#' dle$returnlev
#' 
#' # Estimate discharge that could occur every 80 years (at least empirically):
#' Q80 <- distLextreme(dlf=dle, RPs=80)$returnlev
#' round(sort(Q80[,1], decr=TRUE),1)
#' # 99 to 143 m^3/s can make a relevant difference in engineering!
#' # That's why the rows weighted by GOF are helpful! Weights are given as in
#' distLgofPlot(dle, ranks=FALSE) # See also section weighted mean below
#' # For confidence intervals see distLexBoot
#' 
#' 
#' # Return period of a given discharge value, say 120 m^3/s:
#' sort(1/(1-sapply(dle$parameter, plmomco, x=120) )  )[1:13]
#' # exponential:                 every 29 years
#' # gev (general extreme value dist):  58,
#' # Weibull:                     every 72 years only
#' 
#' 
#' 
#' # Advanced options -------------------------------------------------------------
#' distLextremePlot(dlf=dle)
#' # Line colors / select distributions to be plotted:
#' distLextremePlot(dle, nbest=17, coldist=heat.colors(17), lty=1:5) # lty is recycled
#' distLextremePlot(dle, selection=c("gev", "gam", "gum"), coldist=4:6, PPcol=3, lty=3:2)
#' distLextremePlot(dle, selection=c("gpa","glo","wei","exp"), pch=c(NA,NA,6,8), 
#'                  order=TRUE, cex=c(1,0.6, 1,1), log=TRUE, PPpch=c(16,NA), n_pch=20)
#' # use n_pch to say how many points are drawn per line (important for linear axis) 
#' 
#' # Why do I not get the following warning in interactive mode?
#' ## Warning in if (is.na(dn)) stop("No distributions are left with selection ",  :
#' ## the condition has length > 1 and only the first element will be used
#' # apparently, warnings do not get passed from one function to the next...
#' 
#' distLextremePlot(dle, legarg=list(cex=0.5, x="bottom", box.col="red", col=3))
#' # col in legarg list is (correctly) ignored


#' \dontrun{
#' ## Excluded from package R CMD check because it's time consuming
#' 
#' distLextremePlot(dle, PPpch=c(1,NA)) # only Weibull plotting positions
#' # add different dataset to existing plot:
#' distLextreme(Nile/15, add=TRUE, PPpch=NA, coldist=1, selection="wak", legend=FALSE)
#' 
#' # Logarithmic axis
#' distLextreme(Nile, log=TRUE, nbest=8)
#' 
#' 
#' 
#' # weighted mean based on Goodness of fit (GOF) ---------------------------------
#' # Add discharge weighted average estimate continuously:
#' distLextremePlot(dle, nbest=17, legend=FALSE)
#' abline(h=115.6, v=50)
#' RP <- seq(1, 70, len=100)
#' DischargeEstimate <- distLextreme(dlf=dle, RPs=RP, plot=FALSE)$returnlev
#' lines(RP, DischargeEstimate["weighted2",], lwd=3, col="orange") 
#' 
#' # Or, on log scale:
#' distLextremePlot(dle, nbest=17, legend=FALSE, log=TRUE)
#' abline(h=115.9, v=50)
#' RP <- unique(round(logSpaced(min=1, max=70, n=200, plot=FALSE),2))
#' DischargeEstimate <- distLextreme(dlf=dle, RPs=RP, plot=FALSE)$returnlev
#' lines(RP, DischargeEstimate["weighted2",], lwd=5) 
#' 
#' 
#' # Effect of data proportion used to estimate GOF -------------------------------
#' # Discharge estimated for 50 years return period
#' Goodness <- function(gofProp)
#' { # ExtremeStatistics
#' ES <- distLextreme(annMax, fitargs=list(gofProp=gofProp), plot=FALSE, quiet=TRUE)
#' ES <- cbind(ES$returnlev[1:nrow(ES$gof),], ES$gof)
#' # simple mean:              # plot(sort(GF))
#' av_simple <- mean(ES[,"RP.50"])     # 114.0178  old: 116.8383
#' # weighted average:
#' av_weight <- sum(ES$weight2 * ES[,"RP.50"]) # 115.2842  old: 117.2784
#' # mean of best 3 distribution functions:
#' av_3best <- mean(ES[1:3, "RP.50"]) # 114.7579   old: 118.3722
#' # most functions underestimate the discharge, if we assume that the weibull
#' # PP method correctly calculates the return period of the highest value
#' c(av_simple=av_simple, av_weight=av_weight, av_3best=av_3best)
#' }
#' 
#' Goodness(1)
#' Goodness(0.2)
#' Proportion <- seq(0.05, 1, len=50)
#' GoF <- sapply(Proportion, Goodness) # takes a while to compute
#' plot(Proportion, GoF[3,], type="l", las=1)
#' lines(Proportion, GoF[2,], col=2)
#' lines(Proportion, GoF[1,], col=3)
#' legend("bottomright", col=1:3, lty=1,
#'    legend=c("mean of 3 dists with lowest RMSE", "weighted average", "simple average"))
#' title(main="The proportion of the data included\nfor calculating RMSE does matter!")
#' 
#' 
#' # Minima -----------------------------------------------------------------------
#' 
#' browseURL("http://nrfa.ceh.ac.uk/data/station/meanflow/39072")
#' qfile <- system.file("extdata/discharge39072.csv", package="berryFunctions")
#' Q <- read.table(qfile, skip=19, header=TRUE, sep=",", fill=TRUE)[,1:2]
#' rm(qfile)
#' colnames(Q) <- c("date","discharge")
#' Q$date <- as.Date(Q$date)
#' plot(Q, type="l")
#' Qmax <- tapply(Q$discharge, format(Q$date,"%Y"), max)
#' distLextreme(Qmax, quiet=TRUE)
#' Qmin <- tapply(Q$discharge, format(Q$date,"%Y"), min)
#' dle <- distLextreme(-Qmin, quiet=TRUE, RPs=c(2,5,10,20,50,100,200,500))
#' distLextremePlot(dle, ylim=c(0,-31), yaxs="i", yaxt="n", ylab="Q annual minimum", nbest=14)
#' axis(2, -(0:3*10), 0:3*10, las=1)
#' -dle$returnlev[c(1:14,21), ]
#' # Some distribution functions are an obvious bad choice for this, so I use
#' # weighted 3: Values weighted by GOF of dist only for the best half.
#' # For the Thames in Windsor, we will likely always have > 9 m^3/s streamflow
#' 
#' 
#' # compare extremeStat with other packages: ---------------------------------------
#' library(extRemes)
#' plot(fevd(annMax))
#' par(mfrow=c(1,1))
#' return.level(fevd(annMax, type="GEV")) # "GP", "PP", "Gumbel", "Exponential"
#' distLextreme(dlf=dle, RPs=c(2,20,100))$returnlev["gev",]
#' # differences are small, but noticeable...
#' # if you have time for a more thorough control, please pass me the results!
#' 
#'  
#' # yet another dataset for testing purposes:
#' Dresden_AnnualMax <- c(403, 468, 497, 539, 542, 634, 662, 765, 834, 847, 851, 873,
#' 885, 983, 996, 1020, 1028, 1090, 1096, 1110, 1173, 1180, 1180,
#' 1220, 1270, 1285, 1329, 1360, 1360, 1387, 1401, 1410, 1410, 1456,
#' 1556, 1580, 1610, 1630, 1680, 1734, 1740, 1748, 1780, 1800, 1820,
#' 1896, 1962, 2000, 2010, 2238, 2270, 2860, 4500)
#' distLextreme(Dresden_AnnualMax)
#' } # end dontrun
#' 
#' @param dat Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.
#' @param dlf List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat! DEFAULT: NULL
#' @param selection Selection of distributions. Character vector with types as in \code{\link[lmomco]{lmom2par}}. DEFAULT: NULL
#' @param truncate Proportion truncated for censored quantile, see \code{\link{distLquantile}}. DEFAULT: 0
#' @param RPs ReturnPeriods for which discharge is estimated. DEFAULT: c(2,5,10,20,50)
#' @param progbars Show progress bars for each loop? DEFAULT: TRUE if n>200
#' @param time \code{\link{message}} execution time? DEFAULT: TRUE
#' @param plot Should the return periods and nbest fitted distributions be plotted by a call to \code{\link{distLextremePlot}}? DEFAULT: TRUE
#' @param fitargs List of arguments passed to \code{\link{distLfit}}, like gofProp. DEFAULT: NULL
#' @param quiet Suppress notes and progbars? DEFAULT: FALSE
#' @param \dots Further arguments passed to \code{\link{distLextremePlot}} like order, lty, lwd, ...
#' 
distLextreme <- function(
dat,
dlf=NULL,
selection=NULL,
truncate=0,
RPs=c(2,5,10,20,50),
progbars=length(dlf$dat)>200,
time=TRUE,
plot=TRUE,
fitargs=NULL,
quiet=FALSE,
... )
{
StartTime <- Sys.time()
if(quiet) progbars <- FALSE
if(any(RPs<1.05) & !quiet) on.exit(message("Note in distLextreme: for RPs=1 rather use min(dat)."), add=TRUE)
# fit distributions and calculate goodness of fits -----------------------------
if( is.null(dlf) )  dlf <- do.call(distLfit, owa(list(dat=dat, 
      datname=deparse(substitute(dat)), truncate=truncate, 
      plot=FALSE, selection=selection, time=FALSE, progbars=progbars, quiet=quiet), 
      fitargs, "dat", "datname", "selection"))
# Emptyness check:
if("error" %in% names(dlf)) 
  {
  on.exit(message("Note in distLextreme: There was an error in distLfit, ",
          "thus I'm not plotting anything: ", dlf$error), add=TRUE)
  plot <- FALSE
  }
# Equality check
if(!missing(dat) & !is.null("dlf")) if(any(dlf$dat != dat) & !quiet)
  on.exit(message("Note in distLextreme: 'dat' differs from 'dlf$dat'. 'dat' is ignored."), add=TRUE)
#
# plot -------------------------------------------------------------------------
if(plot) dlf <- distLextremePlot(dlf=dlf, selection=selection, quiet=quiet, ...)
#
# output (discharge) values at return periods ----------------------------------
returnlev <- distLquantile(dlf=dlf, selection=selection, truncate=truncate, 
                           probs=1-1/RPs, empirical=TRUE, quiet=quiet)
# column names:
colnames(returnlev) <- paste0("RP.", RPs)
# Add to output:
dlf$returnlev <- as.data.frame(returnlev)
if(time & !quiet) on.exit(message("distLextreme execution took ", 
  signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds."), add=TRUE)
return(invisible(dlf))
} # end of function
