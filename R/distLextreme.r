# Extreme value statistics for flood risk estimation
# Berry Boessenkool, 2012 (first draft) - 2014 (main update)

# Input: vector with annual discharge maxima
# Ouput: discharge estimates for given return periods, parameter of distributions, quality of fits, plot with linear axis
# See German RclickHandbuch.wordpress.com Kapitel 15

distLextreme <- function(
dat,          # Vector with extreme values e.g. annual discharge maxima. Ignored if dlf is given.
datname,      # Character string for main and output list.
dlf,          # List as returned by \code{\link{distLfit}}, containing the elements \code{dat, parameter, gof}. Overrides dat!
RPs=c(2,5,10,20,50), # ReturnPeriods for which discharge is estimated
speed=TRUE,   # If TRUE, several computationally intensive time consuming distributions are omitted, for the reasons shown in \code{\link[lmomco]{dist.list()}}
ks=TRUE,      # Include ks.test results in dlf$gof? Computing is much faster when FALSE
selection,    # Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}
gofProp=1,    # Upper proportion of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
progbars=TRUE,# Show progress bars for each loop and \code{\link{message}} execution time?
plot=TRUE,    # Should the return periods and nbest fitted distributions be plotted?
add=FALSE,    # If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one
nbest=5,      # Number of distributions plotted, in order of goodness of fit
xlim,         # X-axis limits. DEFAULT: xlim of plotting positions
ylim,         # Y-lim. DEFAULT: from min to extended max
las=1,        # LabelAxisStyle to orient labels, see {code{link{par}}
main=datname, # Title of plot
xlab="Return Period RP  [a]", # X axis label
ylab="Discharge HQ  [m^3/s]", # Y axis label
col="black",  # Plotting point color
pch=c(16,3),  # point characters for plotting positions after Weibull and Grongorton, respectively
cex=1,        # Character EXpansion of plotting points
coldist=rainbow2(nbest), # Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow
lwd=1,        # Line WiDth of distribution lines
legargs=NULL, # list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd
... )         # Further arguments passed to \code{\link{plot}} and \code{\link{lines}}, like lty, lwd, type, pch, ...
{
StartTime <- Sys.time()
# input checks -----------------------------------------------------------------
if(missing(datname)) datname <- deparse(substitute(dat))
if(!missing(dlf)) if(!missing(dat)) if(dlf$dat != dat)
   warning("'dat' differs from 'dlf$dat'. 'dat' is ignored.")
if(!missing(dlf)) datname <- dlf$datname
if(length(pch)==1) pch[2] <-  if(pch[1]==3) 4 else 3
if(!missing(selection))  nbest <- length(selection)
if(length(coldist) != nbest) warning("Length of coldist (",length(coldist),") was not equal to nbest (",nbest,"). Is now recycled.")
coldist <- rep(coldist, length=nbest)
###require(lmomco) # not necessary anymore. it is listed in 'Imports' now...
###require(berryFunctions) # for rsquare, RMSE, logAxis,
#
# fit distributions and calculate goodness of fits -----------------------------
if( missing(dlf) )  dlf <- distLfit(dat=dat, datname=datname, speed=speed, ks=ks,
       selection=selection, gofProp=gofProp, progbars=progbars, plot=FALSE, time=FALSE)
output <- dlf
# objects from the list
dat <- output$dat; parameter <- output$parameter;  gof <- output$gof
# remove NAs, convert to vector:
dat <- as.numeric( dat[!is.na(dat)]  )
#
# Plotting Positions -----------------------------------------------------------
# Calculate PP according to Weibull and Gringorton
# See chapter 15.2 RclickHandbuch.wordpress.com for differences in PP methods
# They're not used for fitting distributions, but for plotting only
m <- sort(dat) # ascendingly sorted extreme values
n <- length(m);  Rank <- 1:n
# RP = Returnperiod = recurrence interval = 1/P_exceedence = 1/(1-P_nonexc.) :
RPw <- 1/(1-  Rank      /(n+1)     )  # Weibull
RPg <- 1/(1- (Rank-0.44)/(n+0.12)  )  # Gringorton (taken from lmom:::evplot.default)
#
dn <- rownames(gof) # distribution names
if(!missing(selection))
   {names(dn) <- dn ;   dn <- dn[selection] ;   names(dn) <- NULL}
#
# plot -------------------------------------------------------------------------
if(plot){
# Calculate plot limits if not given:
if(missing(ylim)) ylim <- c(min(dat), max(dat)+0.1*diff(range(dat)) )
if(missing(xlim)) xlim <- range(RPw, RPg)
# draw discharges over return periods:
plot(1, type="n", las=las, ylim=ylim, xlim=xlim, main=main, ylab=ylab, xlab=xlab, ...) #RPw, m, cex=cex
# range of discharges:
yval <- seq(from=par("usr")[3], to=par("usr")[4], length=300)
# add nbest distributions:
if(nbest > length(dn)) {nbest <- length(dn)}
coldist <- rep(coldist, length=nbest)
for(i in nbest:1) lines(1/(1-plmomco(yval,parameter[[dn[i]]])), yval, col=coldist[i], lwd=lwd, ...)
# Add plotting positions of actual data:
points(RPw, m, pch=pch[1], cex=cex, ...)
points(RPg, m, pch=pch[2], cex=cex, ...)
box()
# write the names of distributions. - legargs: legend arguments
legdef <- list(
  legend=c(if(!is.na(pch[1])) "Weibull plotting positions",
           if(!is.na(pch[2]))"Gringorten plotting positions",dn[1:nbest]),
  pch=c( pch, rep(NA, nbest)),
  lwd=c(NA,NA,rep(lwd,nbest)),
  col=c(col,col,coldist), x="bottomright", cex=0.7, bty="o")
do.call(legend, args=owa(legdef, legargs, c("legend","pch","lwd","col") ))
} # end if plot
# output -----------------------------------------------------------------------
if( require(pbapply,quietly=TRUE) & progbars ) sapply <- pbapply::pbsapply
if(progbars) print("calculating return levels for return periods:")
returnlev <- sapply(dn, function(d) qlmomco(1-1/RPs, parameter[[d]]))  # as.numeric(
# if length(RPs)==1, returnlev is only a vector,
if(is.null(dim(returnlev)))
     returnlev <- as.matrix(returnlev)        # so convert it to a matrix,
else returnlev <- t(returnlev)                # or else transpose it
# column names:
colnames(returnlev) <- paste0("RP.", RPs)
# add goodness of fit
# Add to output:
output$returnlev <- as.data.frame(returnlev)
if(progbars) message("execution took ", signif(difftime(Sys.time(), StartTime, units="s"),2), " seconds.\n")
return(invisible(output))
} # end of function



# old call to fitting. May save a little calculation time sometimes, but is confusing.
if(FALSE)
{
# fit distributions and calculate goodness of fits -----------------------------
if(  missing(dlf) & missing(gofProp)  )
   dlf <- distLfit(dat=dat, datname=datname, speed=speed, selection=selection, gofProp=gofProp,
                   progbars=progbars, plot=FALSE, time=FALSE)
if( !missing(selection)  )
   dlf <- distLfit(dat=dat, datname=datname, speed=speed, selection=selection, gofProp=gofProp,
                   progbars=progbars, plot=FALSE, time=FALSE)
output <- dlf
# save parameter estimation time if dlf is given, but gofProp also:
if(!missing(gofProp)) output <- distLgof(list(dat=output$dat, gofProp=output$gofProp,
                    parameter=output$parameter), progbars=progbars, plot=FALSE)
}