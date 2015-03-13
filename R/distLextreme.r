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
ks=FALSE,     # Include ks.test results in dlf$gof? Computing is much faster when FALSE
selection=NULL, # Selection of distributions, num or char. Can be negative to leave some out if numeric. char as in \code{\link[lmomco]{lmom2par}}
gofProp=1,    # Upper proportion of \code{dat} to compute goodness of fit (dist / ecdf) with. This enables to focus on the dist tail
progbars=TRUE,# Show progress bars for each loop and \code{\link{message}} execution time?
plot=TRUE,    # Should the return periods and nbest fitted distributions be plotted?
add=FALSE,    # If TRUE, plot is not called before adding lines. This lets you add lines highly customized one by one
nbest=5,      # Number of distributions plotted, in order of goodness of fit
xlim=NULL,    # X-axis limits. DEFAULT: xlim of plotting positions
ylim=NULL,    # Y-lim. DEFAULT: from min to extended max
las=1,        # LabelAxisStyle to orient labels, see \code{\link{par}}
main=datname, # Title of plot
xlab="Return Period RP  [a]", # X axis label
ylab="Discharge HQ  [m^3/s]", # Y axis label
col="black",  # Plotting point color
pch=c(16,3),  # point characters for plotting positions after Weibull and Grongorton, respectively
cex=1,        # Character EXpansion of plotting points
coldist=rainbow2(nbest), # Color for each distribution added with \code{\link{lines}}. DEFAULT: rainbow
lwd=1,        # Line WiDth of distribution lines
legend=TRUE,  # Logical. Add a legend?
legargs=NULL, # list of arguments passed to \code{\link{legend}} except for legend, col, pch, lwd
linargs=NULL, # List of arguments passed to \code{\link{lines}} like lty, lwd, type, pch, ...
... )         # Further arguments passed to \code{\link{plot}} like log="x", xaxt="n", ...
{
StartTime <- Sys.time()
# input checks -----------------------------------------------------------------
if(missing(datname)) datname <- deparse(substitute(dat))
if(!missing(dlf)) if(!missing(dat)) if(dlf$dat != dat)
   warning("'dat' differs from 'dlf$dat'. 'dat' is ignored.")
if(!missing(dlf)) datname <- dlf$datname
###require(lmomco) # not necessary anymore. it is listed in 'Imports' now...
###require(berryFunctions) # for rsquare, RMSE, logAxis,   ditto
#
# fit distributions and calculate goodness of fits -----------------------------
if( missing(dlf) )  dlf <- distLfit(dat=dat, datname=datname, speed=speed, ks=ks,
       selection=selection, gofProp=gofProp, progbars=progbars, plot=FALSE, time=FALSE)
output <- dlf
# objects from the list
dat <- output$dat; parameter <- output$parameter;  gof <- output$gof
# remove NAs, convert to vector:
dat <- as.numeric( dat[!is.na(dat)]  )
# plot -------------------------------------------------------------------------
###if(plot) do.call(distLextremePlot, args=list(dlf=dlf, selection=selection, add=add, nbest=nbest,
###    xlim=xlim, ylim=ylim, las=las, main=main, xlab=xlab, ylab=ylab, col=col,
###    pch=pch, cex=cex, coldist=coldist, lwd=lwd, legend=legend, legargs=legargs, linargs=linargs, ...))
if(plot) output <- distLextremePlot(dlf=dlf, selection=selection, add=add, nbest=nbest,   # output <-
    xlim=xlim, ylim=ylim, las=las, main=main, xlab=xlab, ylab=ylab, col=col,
    pch=pch, cex=cex, coldist=coldist, lwd=lwd, legend=legend, legargs=legargs, linargs=linargs, ...)
# output (discharge) values at return periods ----------------------------------
dn <- rownames(gof) # distribution names
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
