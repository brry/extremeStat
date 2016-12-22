#' Distribution rank comparison
#' 
#' Plot rank comparison of fitted distributions calculated by \code{\link{distLgof}}.
#' 
#' @return None.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{distLgof}}, \code{\link{distLfit}}
#' @keywords hplot distribution
#' @export
#' @importFrom berryFunctions owa
#' @examples
#' # see distLgof
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the element \code{gof}
#' @param ranks plot ranks by different measures of goodness of fit? DEFAULT: TRUE
#' @param weights plot weights in second panel? DEFAULT: TRUE
#' @param nbest Number of distributions plotted by \code{\link{distLplot}}. DEFAULT: NA
#' @param order If TRUE, the best distribution is at the top. False to compare across different dlfs. DEFAULT: TRUE
#' @param add If TRUE, just the lines and legend are added to an existing plot. DEFAULT: FALSE
#' @param type type of line-point combination as in \code{\link{plot}}. DEFAULT: "o"
#' @param col Vector with 3 values for line customization. Recycled if necessary. DEFAULT: c(1,3,4)
#' @param pch Vector with 3 values for line customization. Recycled if necessary. DEFAULT: 1:3
#' @param lty Vector with 3 values for line customization. Recycled if necessary. DEFAULT: 1
#' @param lwd Vector with 3 values for line customization. Recycled if necessary. DEFAULT: 1
#' @param legargs List of arguments passed to \code{\link{legend}} if weights=TRUE, like cex, bg, etc.
#' @param quiet Suppress notes? DEFAULT: FALSE
#' @param main plot titles. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{plot}}.
#' 
distLgofPlot <- function(
dlf,
ranks=TRUE,
weights=TRUE,
nbest=NA,
order=TRUE,
add=FALSE,
type="o",
col=c(1:3,8,4),
pch=c(1:4,NA),
lty=1,
lwd=1,
legargs=NULL,
quiet=FALSE,
main=NULL,
...)
{
# Input control:
if(nrow(dlf$gof)<1) stop("No fitted distributions in dlf.")
if(nrow(dlf$gof)<2) stop("Only ",toString(rownames(dlf$gof)), " was fitted, thus GOF can't be compared.")
# recycling:
col <- rep(col, length=5)
pch <- rep(pch, length=5)
lty <- rep(lty, length=5)
lwd <- rep(lwd, length=5)
# Objects from list:
gof <- dlf$gof
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
if(is.null(main2)) main2 <- paste("Ranking of distributions by goodness of fit")
#
if(!add) plot(Ranks[,1], 1:n, type="n", ylab="", xlab="Rank", yaxt="n", main=main2, ...)
if(!add) axis(2, 1:n, rownames(Ranks), las=2)
lines(Ranks[,1], 1:n, pch=pch[1], col=col[1], lty=lty[1], lwd=lwd[1], type=type)
lines(Ranks[,2], 1:n, pch=pch[2], col=col[2], lty=lty[2], lwd=lwd[2], type=type)
if(ks)
lines(Ranks[,4], 1:n, pch=pch[3], col=col[3], lty=lty[3], lwd=lwd[3], type=type)
abline(h=n-nbest+0.5)
legend("bottomleft", c("RMSE","R2",if(ks)"ks.test"), lty=lty, col=col, pch=pch)
}
# Weights:  -------------------------------------------------------------------
if(weights)
{
if(is.null(main)) main <- "Weights of distributions\nfor weighted average"
Xlim <- range(gof[, grep("weight", colnames(gof))], na.rm=T)
plot(gof$weightc, nrow(gof):1, type="n", xlim=Xlim, yaxt="n", xlab="Weight",
     las=1, ylab="Dist", main=main)
lines(gof$weight1, nrow(gof):1, pch=pch[1], col=col[1], lty=lty[1], lwd=lwd[1], type=type)
lines(gof$weight2, nrow(gof):1, pch=pch[2], col=col[2], lty=lty[2], lwd=lwd[2], type=type)
lines(gof$weight3, nrow(gof):1, pch=pch[3], col=col[3], lty=lty[3], lwd=lwd[3], type=type)
lines(gof$weightc, nrow(gof):1, pch=pch[4], col=col[4], lty=lty[4], lwd=lwd[4], type=type)
lines(gof$RMSE,    nrow(gof):1, pch=pch[5], col=col[5], lty=lty[5], lwd=lwd[5], type=type)
text(head(gof$RMSE,1), nrow(gof), "RMSE", col=col[5], adj=-0.2)
axis(2, nrow(gof):1, rownames(gof), las=2)
do.call(legend, berryFunctions::owa(list(x="bottomright", legend=c("1: max(r)-r + min(r)",
        "2: max(r)-r", "3: first half", "c: custom"),
       title="Weighted by RMSE = r", col=col[1:4], pch=pch[1:4], lty=lty[1:4], lwd=lwd[1:4]), legargs))
}
} # end of function
