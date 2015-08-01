# Extreme value statistics for flood risk estimation
# Berry Boessenkool, 2015

distLextremePlot <- function(
dlf,          # List as returned by \code{\link{distLextreme}}, containing the elements \code{dat, parameter, gof}. Overrides dat!
selection=NULL,# Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}
order=FALSE,  # If selection is given, should legend and colors be ordered by gof anyways?
add=FALSE,    # If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one
nbest=5,      # Number of distributions plotted, in order of goodness of fit
xlim=NULL,         # X-axis limits. DEFAULT: xlim of plotting positions
ylim=NULL,         # Y-lim. DEFAULT: from min to extended max
las=1,        # LabelAxisStyle to orient labels, see \code{\link{par}}
main=dlf$datname, # Title of plot
xlab="Return Period RP  [a]", # X axis label
ylab="Discharge HQ  [m^3/s]", # Y axis label
col="black",  # Plotting point colors, vector of length two for Weibull and Gringorton, recycled
pch=c(16,3),  # point characters for plotting positions after Weibull and Grongorton, respectively
cex=1,        # Character EXpansion of plotting points
coldist=rainbow2(nbest), # Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow
lty=1,        # Line TYpe for plotted distributions. Recycled vector of length nbest.  
lwd=1,        # Line WiDth of distribution lines. Not recycled.
legend=TRUE,  # Logical. Add a legend?
legargs=NULL, # list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd
linargs=NULL, # List of arguments passed to \code{\link{lines}} like lty, lwd, type, pch, ...
quiet=FALSE,  # Suppress notes?
... )         # Further arguments passed to \code{\link{plot}} like log="x", xaxt="n", ...
{
# Input operations:
output <- dlf
if(length(pch)==1) pch[2] <- if(is.na(pch)) NA else if(pch[1]==3) 4 else 3
col <- rep(col, length.out=2)
if(!is.null(selection))  nbest <- length(selection)
if(length(coldist) != nbest & !quiet)
   message("note in distLextremePlot: Length of coldist (",length(coldist),") was not equal to nbest (",nbest,"). Is now recycled.")
coldist <- rep(coldist, length=nbest)
# Extract objects from dlf:
dat <- dlf$dat
parameter <- dlf$parameter
gof <- dlf$gof
# Plotting Positions -----------------------------------------------------------
# Calculate PP according to Weibull and Gringorton
# See chapter 15.2 RclickHandbuch.wordpress.com for differences in PP methods
# They're not used for fitting distributions, but for plotting only
m <- sort(dat) # ascendingly sorted extreme values
n <- length(m);  Rank <- 1:n
# RP = Returnperiod = recurrence interval = 1/P_exceedence = 1/(1-P_nonexc.) :
RPw <- 1/(1-  Rank      /(n+1)     )  # Weibull
RPg <- 1/(1- (Rank-0.44)/(n+0.12)  )  # Gringorton (taken from lmom:::evplot.default)
#  Selection -------------------------------------------------------------------
dn <- rownames(gof) # distribution names
if(!is.null(selection))
  {
  names(dn) <- dn
  selection <- selection[selection %in% dn]
  dn <- dn[selection]
  if(order)
    {
    dno <- rownames(gof)
    dno <- dno[dno %in% dn]
    dn <- dn[dno]
    }
  names(dn) <- NULL
  if(all(is.na(dn))) stop("No distributions are left with selection ", paste(selection, collapse=", "))
  }
# Plotting ---------------------------------------------------------------------
# Calculate plot limits if not given:
if(is.null(ylim)) ylim <- c(min(dat), max(dat)+0.1*diff(range(dat)) )
if(is.null(xlim)) xlim <- range(RPw, RPg)
# draw discharges over return periods:
if(!add) plot(1, type="n", las=las, ylim=ylim, xlim=xlim, main=main, ylab=ylab, xlab=xlab, ...) #RPw, m, cex=cex
# range of discharges:
yval <- seq(from=par("usr")[3], to=par("usr")[4], length=300)
# add nbest distributions as lines:
if(nbest > length(dn)) {nbest <- length(dn)}
coldist <- rep(coldist, length=nbest)
lty <- rep(lty, length=nbest)
for(i in nbest:1)
  {
  Pnonexceed <- plmomco(yval,parameter[[dn[i]]]) # print(Pnonexceed, digits=20)
  Pnonexceed[Pnonexceed>1] <- 1 # remove numerical errors
  linargsdefault <- list(x=1/(1-Pnonexceed), y=yval, col=coldist[i], lty=lty[i], lwd=lwd)
  do.call(lines, args=owa(linargsdefault, linargs))
  }
# Add plotting positions of actual data:
points(RPw, m, pch=pch[1], cex=cex, col=col[1])
points(RPg, m, pch=pch[2], cex=cex, col=col[2])
box()
# Legend -----------------------------------------------------------------------
# write the names of distributions. - legargs: legend arguments
if(legend){
legdef <- list(
  legend=c(if(!is.na(pch[1])) "Weibull plotting positions",
           if(!is.na(pch[2]))"Gringorten plotting positions",dn[1:nbest]),
  pch=c( pch, rep(NA, nbest)),
  lwd=c(if(!is.na(pch[1])) NA,     if(!is.na(pch[2])) NA,     rep(lwd,nbest)),
  col=c(if(!is.na(pch[1])) col[1], if(!is.na(pch[2])) col[2], coldist),
  lty=c(if(!is.na(pch[1])) NA,     if(!is.na(pch[2])) NA,     lty),
  x="bottomright", cex=0.7, bty="o")
do.call(graphics::legend, args=owa(legdef, legargs, "legend","pch","lwd","col","lty"))
}
# output dlf object
output$RPweibull <- RPw
output$RPgringorton <- RPg
return(invisible(output))
}
