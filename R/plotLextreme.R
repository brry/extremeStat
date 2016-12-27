#' Plot extreme value statistics
#' 
#' Plots distributions fitted by L-moments and adds plotting positions by Weibull and Gringorton.
#'
#' @details This is an auxiliary graphing function to \code{\link{distLextreme}}
#' 
#' @return none, plots things.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, March 2015, updated heavily Aug 2015
#' @seealso \code{\link{distLextreme}}, \code{\link{plotLfit}}
#' @keywords hplot dplot distribution
#' @export
#' @importFrom berryFunctions rainbow2 owa logAxis
#' @importFrom lmomco plmomco
#' 
#' @examples
#' #see
#' ?distLextreme
#' 
#' @param dlf     List as returned by \code{\link{distLextreme}}, containing the 
#'                elements \code{dat, dleB <- distLexBoot(dlf, nbest=4, conf.lev=0.5), gof}.
#' @param selection Selection of distributions. Character vector with type as in 
#'                \code{\link[lmomco]{lmom2par}}. DEFAULT: NULL
#' @param order   If selection is given, should legend and colors be ordered 
#'                by gof anyways? DEFAULT: FALSE
#' @param add     If TRUE, plot is not called before adding lines. This lets you 
#'                add lines to an existing plot. DEFAULT: FALSE
#' @param nbest   Number of distributions plotted, in order of goodness of fit. 
#'                Overwritten internally if selection is given. DEFAULT: 5
#' @param log     Logarithmic x-axis? DEFAULT: FALSE
#' @param xlim    X-axis limits. DEFAULT: xlim of plotting positions
#' @param ylim    Y-lim. DEFAULT: from min to extended max
#' @param las     LabelAxisStyle to orient labels, see \code{\link{par}}. DEFAULT: 1
#' @param main    Title of plot. DEFAULT: dlf$datname
#' @param xlab    X axis label. DEFAULT: "Return Period RP  [a]"
#' @param ylab    Y axis label. Please note that the ubuntu pdf viewer might be 
#'                unable to display unicode superscript. DEFAULT: "Discharge HQ  [m3/s]"
#' @param PPcol   Plotting Position point colors, vector of length two for 
#'                Weibull and Gringorton, recycled. PP are not used for 
#'                fitting distributions, but for plotting only. DEFAULT: "black"
#' @param PPpch   point characters for plotting positions after Weibull and 
#'                Gringorton, respectively. NA to suppress in plot and legend. DEFAULT: c(16,3)
#' @param PPcex   Character EXpansion of plotting points. DEFAULT: 1
#' @param coldist Color for each distribution added with \code{\link{lines}}. 
#'                Recycled, if necessary. DEFAULT: \code{\link[berryFunctions]{rainbow2}}
#' @param lty     Line TYpe for plotted distributions. Is recycled to from a 
#'                vector of length nbest, i.e. a value for each dist. DEFAULT: 1
#' @param lwd     Line WiDth of distribution lines. Recycled vector of length nbest. 
#'                DEFAULT: 1
#' @param pch     Point CHaracter of points added at regular intervals. 
#'                NA to suppress. Recycled vector of length nbest. DEFAULT: NA
#' @param cex     if pch != NA, size of points. Recycled vector of length nbest. 
#'                DEFAULT: 1
#' @param n_pch   Number of points spread evenly along the line. 
#'                Recycled vector of length nbest. DEFAULT: 15
#' @param legend  Logical. Add a legend? DEFAULT: TRUE
#' @param legargs list of arguments passed to \code{\link{legend}} except for 
#'                legend, col, pch, lwd, lty. DEFAULT: NULL
#' @param quiet   Suppress notes? DEFAULT: FALSE
#' @param \dots   Further arguments passed to \code{\link{plot}} like yaxt="n", ...
#' 
plotLextreme <- function(
dlf,
selection=NULL,
order=FALSE,
add=FALSE,
nbest=5,
log=FALSE,
xlim=NULL,
ylim=NULL,
las=1,
main=dlf$datname,
xlab="Return Period RP  [a]",
ylab="Discharge HQ  [m\U00B3/s]",
PPcol="black",
PPpch=c(16,3),
PPcex=1,
coldist=berryFunctions::rainbow2(nbest),
lty=1,
lwd=1,
pch=NA,
cex=1,
n_pch=15,
legend=TRUE,
legargs=NULL,
quiet=FALSE,
... )
{
# PP (Plotting Position) points charactereistics recycled:
if(length(PPpch)==1) PPpch[2] <- if(is.na(PPpch)) NA else if(PPpch[1]==3) 4 else 3
PPcol <- rep(PPcol, length.out=2)
PPcex <- rep(PPcex, length.out=2)
# Plotting Positions -----------------------------------------------------------
# Calculate PP according to Weibull and Gringorton
# See chapter 15.2 RclickHandbuch.wordpress.com for differences in PP methods
# They're not used for fitting distributions, but for plotting only
m <- sort(dlf$dat) # ascendingly sorted extreme values
n <- length(m);  Rank <- 1:n
# RP = Returnperiod = recurrence interval = 1/P_exceedence = 1/(1-P_nonexc.) :
RPw <- 1/(1-  Rank      /(n+1)     )  # Weibull
RPg <- 1/(1- (Rank-0.44)/(n+0.12)  )  # Gringorton (taken from lmom:::evplot.default)
#  Selection -------------------------------------------------------------------
dn <- rownames(dlf$gof)[order(dlf$gof$RMSE)] # distribution names
if(!is.null(selection))
  {
  names(dn) <- dn
  sind <- selection %in% dn
  if(!any(sind)) stop("selection ", toString(selection), 
                                   " is not available in dlf$gof.")
  if(any(!sind)) message("Note in plotLextreme: selection ", toString(selection[!sind]), 
                         " is not available in dlf$gof, thus ignored.")
  selection <- selection[sind]
  dn <- dn[selection]
  if(order)
    {
    dno <- rownames(dlf$gof)
    dno <- dno[dno %in% dn]
    dn <- dn[dno] # in descending order of goodness of fit
    }
  names(dn) <- NULL
  if(all(is.na(dn))) stop("No distributions are left with current selection.")
  # nbest:
  nbest <- length(dn) # shorten to selection
  }
# nbest must be shorter if dlf with fewer parameters is given as input:
nbest <- pmin(length(dlf$parameter), nbest)
# control coldist length
if(length(coldist) != nbest & !quiet)
  {
  # This happens of coldist is specified with wrong length.
  # Can happen if selection is truncated (misspellings, dists not fitted)
  message("Note in plotLextreme: Length of coldist (",length(coldist),
          ") was not equal to nbest (",nbest,"). Is now recycled.")
  coldist <- rep(coldist, length=nbest)
  }
#
# Plotting ---------------------------------------------------------------------
# Calculate plot limits if not given:
if(is.null(ylim)) ylim <- c(min(dlf$dat), max(dlf$dat)+0.1*diff(range(dlf$dat)) )
if(is.null(xlim)) xlim <- range(RPw, RPg)
# draw discharges over return periods:
if(!add) plot(1, type="n", las=las, ylim=ylim, xlim=xlim, main=main, ylab=ylab, xlab=xlab, 
              log=if(log) "x" else "", xaxt=if(log) "n" else "s", ...)
if(log) berryFunctions::logAxis(1)
# range of discharges:
yval <- seq(from=par("usr")[3], to=par("usr")[4], length=300)
y_int <- approx(yval, n=n_pch)$y
# Rycycle distfun lines arguments:
lty <- rep(lty, length=nbest)
lwd <- rep(lwd, length=nbest)
pch <- rep(pch, length=nbest)
cex <- rep(cex, length=nbest)
# add nbest distributions as lines:
for(i in nbest:1)
  {
  Pnonexceed <- lmomco::plmomco(yval,dlf$parameter[[dn[i]]]) # print(Pnonexceed, digits=20)
  Pnonexceed[Pnonexceed>1] <- 1 # remove numerical errors
  xval <- 1/(1-Pnonexceed)
  lines(x=xval, y=yval, col=coldist[i], lty=lty[i], lwd=lwd[i])
  if(!is.na(pch[i]))
    {
    x_int <- approx(xval, n=n_pch)$y
    points(x_int, y_int, pch=pch[i], col=coldist[i])
    }
  }
# Add plotting positions of actual data:
points(RPw, m, pch=PPpch[1], cex=PPcex, col=PPcol[1])
points(RPg, m, pch=PPpch[2], cex=PPcex, col=PPcol[2])
box()
# Legend -----------------------------------------------------------------------
# write the names of distributions. - legargs: legend arguments
if(legend){
PP1 <- !is.na(PPpch[1])
PP2 <- !is.na(PPpch[2])
legdef <- list(
  legend=c(if(PP1) "Weibull PP",if(PP2)"Gringorten PP", dn[1:nbest]),
  pch=   c(if(PP1) PPpch[1],    if(PP2) PPpch[2],       pch),
  pt.cex=c(if(PP1) 1,           if(PP2) 1,              cex),
  lwd=   c(if(PP1) NA,          if(PP2) NA,             lwd),
  col=   c(if(PP1) PPcol[1],    if(PP2) PPcol[2],       coldist),
  lty=   c(if(PP1) NA,          if(PP2) NA,             lty),
  x="bottomright",  
  cex=0.8, bg="white")
do.call(graphics::legend, args=berryFunctions::owa(legdef, legargs, 
                                                   "legend","pch","lwd","col","lty"))
}
# output dlf object
dlf$RPweibull <- RPw
dlf$RPgringorton <- RPg
dlf$coldist <- coldist
return(invisible(dlf))
}
