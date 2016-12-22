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
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' # see distLgof
#' 
#' @param dlf List as returned by \code{\link{distLfit}}, containing the element \code{gof}
#' @param type,col,pch,lty,lwd Vectors with 5 values for line customization. Recycled if necessary.
#' @param legargs List of arguments passed to \code{\link{legend}}, like cex, bg, etc.
#' @param main,xlab,ylab plot title and axis labels
#' @param xlim Range of x axis. DEFAULT: range(gof$weight*)
#' @param \dots Further arguments passed to \code{\link{plot}}.
#' 
plotLgof <- function(
dlf,
type="o",
col=RColorBrewer::brewer.pal(5,"Set2"),#c("#66C2A5", "#FC8D62", "#E78AC3",8,4),
pch=c(1:4,NA),
lty=1,
lwd=1,
legargs=NULL,
main="Distribution function GOF and weights",
xlab="Weight",
ylab="",
xlim=range(gof[,grep("weight",colnames(gof))], na.rm=TRUE),
...)
{
# Object from list, for code readability:
gof <- dlf$gof

# Input control:
if(nrow(gof)<1) stop("No fitted distributions in dlf.")
if(nrow(gof)<2) stop("Only ", toString(rownames(gof)), 
                     " was fitted, thus GOF can't be compared.")
# recycling:
pch <- rep(pch, length=5)
col <- rep(col, length=5)
lty <- rep(lty, length=5)
lwd <- rep(lwd, length=5)
type<- rep(type,length=5)

# plotting
plot(1, type="n", xlim=xlim, ylim=c(1,nrow(gof)), yaxt="n", xlab=xlab, ylab=ylab, main=main, ...)
cnames <- c("weight1","weight2","weight3","weightc","RMSE")
for(i in 1:5) lines(gof[,cnames[i]], nrow(gof):1, 
                    pch=pch[i], col=col[i], lty=lty[i], lwd=lwd[i], type=type[i])
text(gof$RMSE[1], nrow(gof), "RMSE", col=col[5], adj=-0.2)
axis(2, nrow(gof):1, rownames(gof), las=1)
do.call(legend, berryFunctions::owa(list(x="bottomright", 
        legend=c("1: max(r)-r + min(r)", "2: max(r)-r", "3: first half", "c: custom"),
        title="Weighted by RMSE = r", 
        pch=pch[1:4], col=col[1:4], lty=lty[1:4], lwd=lwd[1:4]), legargs))

} # end of function
