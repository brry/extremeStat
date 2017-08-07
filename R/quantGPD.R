#' @title Fast GPD quantile estimate
#' @description Fast GPD quantile estimate through L-moments
#' @return Vector with quantiles
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jun 2017
#' @seealso \code{\link{q_gpd}} for a comparison across R packages and methods, \code{\link{distLquantile}} to compare distributions
#' @keywords distribution robust univar
#' @importFrom lmomco lmoms pargpa qlmomco
#' @export
#' @examples
#' data(annMax)
#' quantile(annMax, 0.99)
#' quantGPD(annMax, 0.99)
#' 
#' \dontrun{ # Excluded from CRAN checks to reduce checking time
#' data(rain, package="ismev") ;  rain <- rain[rain>0]
#' hist(rain, breaks=50, col=7)
#' tr <- seq(0,0.999, len=50)
#' qu <- pbapply::pbsapply(tr, quantGPD, x=rain, probs=c(0.9,0.99,0.999) ) # 30 s
#' plot(tr, qu[3,], ylim=range(rain), las=1, type="l")
#' lines(tr, qu[2,], col=2); lines(tr, qu[1,], col=4)
#' 
#' tr <- seq(0.88,0.999, len=50)
#' qu <- pbapply::pbsapply(tr, quantGPD, x=rain, probs=c(0.9,0.99,0.999) ) # 5 s
#' plot(tr, qu[3,], ylim=range(rain), las=1, type="l")
#' lines(tr, qu[2,], col=2); lines(tr, qu[1,], col=4);
#' tail(qu["n",])
#' 
#' library(microbenchmark)
#' data(rain, package="ismev"); rain <- rain[rain>0]
#' mb <- microbenchmark(quantGPD(rain[1:200], truncate=0.8, probs=0.99, addn=F),
#' distLquantile(rain[1:200], sel="gpa", emp=F, truncate=0.8, quiet=T, probs=0.99)[1,1]
#' )
#' boxplot(mb)
#' # since computing the lmoments takes most of the computational time,
#' # there's not much to optimize in large samples like n=2000
#' 
#' }
#' 
#' @param x         Vector with numeric values. NAs are silently ignored.
#' @param probs     Probabilities. DEFAULT: c(0.8,0.9,0.99)
#' @param truncate,threshold Truncation proportion or threshold. DEFAULT: 0, computed
#'                  See \code{\link{q_gpd}}.
#' @param quiet     Should messages from this function be suppressed? DEFAULT: FALSE
#' @param addn      Logical: add element with sample size (after truncation). DEFAULT: TRUE
#' @param \dots     Further arguments passed to \code{lmomco::\link[lmomco]{pargpa}}
#' 
quantGPD <- function(
x,
probs=c(0.8,0.9,0.99),
truncate=0,
threshold=berryFunctions::quantileMean(x, truncate),
addn=TRUE,
quiet=FALSE,
...
)
{
x <- as.numeric(x)
x <- x[is.finite(x)]
x <- x[x>=threshold]
mom <- lmomco::lmoms(x, nmom=5)
par <- lmomco::pargpa(mom, ...)
probs2 <- probs
if(truncate!=0)
  {
  probs2 <- (probs-truncate)/(1-truncate)
  probs2[probs < truncate] <- 0
  }
out <- tryStack(lmomco::qlmomco(f=probs2, para=par))
names(out) <- paste0(probs*100,"%")
out[probs < truncate] <- NA
if(addn) out <- c(out, n=length(x))
out
}
