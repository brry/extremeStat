# distributions via linear moments
# Berry Boessenkool, Sept 2014
# compare ranks of 'Goodness of Fit' measures in distLfit

#' rank comparison
#' 
#' plot rank comparison of fitted distributions calculated by \code{\link{distLgof}}.
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the elements \code{gof, gofProp}
#' @param ranks plot ranks by different measures of goodness of fit? DEFAULT: TRUE
#' @param weights plot weights in second panel? DEFAULT: TRUE
#' @param nbest Number of distributions plotted by \code{\link{distLplot}}. DEFAULT: NA
#' @param order If TRUE, the best distribution is at the top. False to compare across different dlfs. DEFAULT: TRUE
#' @param add If TRUE, just the lines and legend are added to an existing plot. DEFAULT: FALSE
#' @param type type of line-point combination as in \code{\link{plot}}. DEFAULT: "o"
#' @param col Vector with 3 values for line customization. Recycled if necessary. DEFAULT: c(1,3,4)
#' @param pch Vector with 3 values for line customization. Recycled if necessary. DEFAULT: 1:3
#' @param lty Vector with 3 values for line customization. Recycled if necessary. DEFAULT: 1
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param main plot titles. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{plot}}.

#' @return None.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{distLgof}}, \code{\link{distLfit}}
#' @keywords hplot distribution
#' @export
#' @importFrom berryFunctions pastec
#' @examples
#' # see distLgof
#' 
distLgofPlot <- function(
dlf, # List as returned by \code{\link{distLfit}}, containing the elements \code{gof, gofProp}
ranks=TRUE,
weights=TRUE,
nbest=NA,
order=TRUE,
add=FALSE,
type="o",
col=c(1,3,4),
pch=1:3,
lty=1,
quiet=FALSE,
main=NULL,
...)
{
# Input control:
if(nrow(dlf$gof)<1) stop("No fitted distributions in dlf.")
if(nrow(dlf$gof)<2) stop("Only ", berryFunctions::pastec(rownames(dlf$gof)), " was fitted, thus GOF can't be compared.")
# recycling:
col <- rep(col, length=3)
pch <- rep(pch, length=3)
lty <- rep(lty, length=3)
# Objects from list:
gof <- dlf$gof
gofProp <- dlf$gofProp
# catch ks=FALSE results:
ks <- "ksP" %in% colnames(gof)
if(!ks & !quiet) on.exit(message(
  "note in distLgofPlot: This result from distLgof with ks=FALSE ignores ks ranks."), add=TRUE)
if(ranks & weights) { op <- par(mfrow=c(1,2)) ;  on.exit(par(op), add=TRUE) }
# Ranks: -----------------------------------------------------------------------
Ranks <- cbind(rank(gof$RMSE), rank(-gof$R2))
if(ks) Ranks <- cbind(Ranks, rank(-gof$ksP), rank(gof$ksD))
rownames(Ranks) <- rownames(gof)
if(order) Ranks <- Ranks[ order(Ranks[,1],     decreasing=TRUE) , ]  #order(rowMeans(Ranks)
else      Ranks <- Ranks[ order(rownames(gof), decreasing=TRUE) , ]
n <- nrow(Ranks)
# Plot:
if(ranks)
{
main2 <- main
if(is.null(main2)) main2 <- paste("Ranking of distributions by goodness of fit\n",
                               "Proportion of data for RMSE/R2 :", gofProp*100, "%")
#
if(!add) plot(Ranks[,1], 1:n, type="n", ylab="", xlab="Rank", yaxt="n", main=main2, ...)
if(!add) axis(2, 1:n, rownames(Ranks), las=2)
lines(Ranks[,1], 1:n, type=type, col=col[1], pch=pch[1], lty=lty[1])
lines(Ranks[,2], 1:n, type=type, col=col[2], pch=pch[2], lty=lty[2])
if(ks)
lines(Ranks[,4], 1:n, type=type, col=col[3], pch=pch[3], lty=lty[3])
abline(h=n-nbest+0.5)
legend("bottomleft", c("RMSE","R2",if(ks)"ks.test"), lty=lty, col=col, pch=pch)
}
# Weights:  -------------------------------------------------------------------
if(weights)
{
if(is.null(main)) main <- "Weights of distributions\nfor weighted average"
Xlim <- range(gof[, grep("weight", colnames(gof))], na.rm=T)
plot(gof$weightc, nrow(gof):1, xlim=Xlim, type="o", yaxt="n", xlab="Weight",
     las=1, ylab="Dist", main=main, pch=4, col=8)
lines(gof$weight1, nrow(gof):1, pch=1, col=1, type="o")
lines(gof$weight2, nrow(gof):1, pch=2, col=2, type="o")
lines(gof$weight3, nrow(gof):1, pch=3, col=3, type="o")
lines(gof$RMSE,    nrow(gof):1, col=4)
text(head(gof$RMSE,1), nrow(gof), "RMSE", col=4, adj=-0.2)
axis(2, nrow(gof):1, rownames(gof), las=2)
legend("bottomright", c("1: max(r)-r + min(r)", "2: max(r)-r", "3: first half", "c: custom"),
       title="Weighted by RMSE = r", col=c(1:3,8), pch=1:4, lty=1)
}
} # end of function
