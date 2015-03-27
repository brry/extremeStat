# distributions via linear moments
# Berry Boessenkool, Sept 2014
# compare ranks of 'Goodness of Fit' measures in distLfit

distLgofPlot <- function(
dlf, # List as returned by \code{\link{distLfit}}, containing the elements \code{gof, gofProp}
nbest=NA, # Number of distributions plotted by \code{\link{distLplot}}
order=TRUE, # If TRUE, the best distribution is at the top. False to compare across different dlfs.
add=FALSE, # If TRUE, just the lines and legend are added to an existing plot.
type="o", # type of line-point combination as in \code{\link{plot}}
col=c(1,3,4), pch=1:3, lty=1, # Vectors with 3 values for line customization. Recycled if necessary. DEFAULTS: see created plot
quiet=FALSE, # Suppress notes?
...) # Further arguments passed to \code{\link{plot}}.
{
# recycling:
col <- rep(col, length=3)
pch <- rep(pch, length=3)
lty <- rep(lty, length=3)
# Objects from list:
gof <- dlf$gof
gofProp <- dlf$gofProp
# catch ks=FALSE results:
ks <- "ksP" %in% colnames(gof)
if(!ks & !quiet) message("note in distLgofPlot: This result from distLgof with ks=FALSE ignores ks ranks.")
# Ranks:
Ranks <- cbind(rank(gof$RMSE), rank(-gof$R2))
if(ks) Ranks <- cbind(Ranks, rank(-gof$ksP), rank(gof$ksD))
rownames(Ranks) <- rownames(gof)
if(order) Ranks <- Ranks[ order(Ranks[,1],     decreasing=TRUE) , ]  #order(rowMeans(Ranks)
else      Ranks <- Ranks[ order(rownames(gof), decreasing=TRUE) , ]
n <- nrow(Ranks)
# Plot:
if(!exists("main")) main <- paste("Ranking of distributions by goodness of fit\n",
                               "Proportion of data for RMSE/R2 :", gofProp*100, "%")
#
if(!add) plot(Ranks[,1], 1:n, type="n", ylab="", xlab="Rank", yaxt="n", main=main, ...)
if(!add) axis(2, 1:n, rownames(Ranks), las=2)
lines(Ranks[,1], 1:n, type=type, col=col[1], pch=pch[1], lty=lty[1])
lines(Ranks[,2], 1:n, type=type, col=col[2], pch=pch[2], lty=lty[2])
if(ks)
lines(Ranks[,4], 1:n, type=type, col=col[3], pch=pch[3], lty=lty[3])
abline(h=n-nbest+0.5)
legend("bottomleft", c("RMSE","R2",if(ks)"ks.test"), lty=lty, col=col, pch=pch)
} # end of function
