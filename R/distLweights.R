#' Compute distribution weights from GOF
#'
#' Determine distribution function weights from RMSE for weighted averages.
#' The weights are inverse to RMSE, weight1 for all dists, 
#' weight2 places zero weight on the worst function, weight3 on the worst half of functions.
#'
#' @return data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLgof}}, \code{\link{distLquantile}}
#' @keywords distribution
#' @export
#' @examples
#' 
#' distLweights(c(gum=0.20, wak=0.17, gam=0.21))
#' distLweights(c(gum=0.20, wak=0.17, gam=0.21), order=FALSE)
#' distLweights(c(gum=0.20))
#' df <- data.frame(gum=2:5, rmse=3:6)
#' rownames(df) <- letters[1:4]
#' distLweights(df)
#' 
#' set.seed(42); x <- data.frame(A=1:5, RMSE=runif(5)) ; x
#' distLweights(x)
#' distLweights(x, weightc=c("1"=3, "3"=5)) 
#' distLweights(x, weightc=c("1"=3, "3"=5), order=FALSE) 
#' 
#' distLweights(data.frame(rmse=1:2))
#' distLweights(data.frame(Rmse=1:3))
#' distLweights(data.frame(rmse=1:6))
#' distLweights(data.frame(rmse=1:11))
#' distLweights(data.frame(rmse=1:12))
#' distLweights(data.frame(rmse=sample(1:12)))
#' distLweights(data.frame(rmse=sample(1:12)), order=FALSE)
#' 
#' 
#' # RMSE vs R2 for GOF judgement --------
#' library(lmomco)
#' data(annMax)
#' 
#' compranks <- function(d)
#' {
#' x <- sort(annMax, decreasing=TRUE)
#' tcdfs <- plmomco(x,dlf$parameter[[d]])
#' ecdfs <- ecdf(annMax)(x) # Empirical CDF
#' # Root Mean Square Error, R squared:
#' berryFunctions::linReg(tcdfs, ecdfs, lwd=1, pch=16, main=d, digits=5, xlim=c(0.5, 1),
#'        ylim=c(0.5, 1), pos1="topleft")
#' abline(a=0,b=1, lty=3)
#' c(berryFunctions::rmse(tcdfs, ecdfs), berryFunctions::rsquare(tcdfs, ecdfs))
#' }
#' 
#' dlf <- distLfit(annMax)
#' dn <- rownames(dlf$gof)
#' 
#' op <- par(mfrow=c(5,4), mar=rep(1,4), xaxt="n", yaxt="n")
#' for(i in dn) compranks(i)
#' par(op)
#' # so revgum, nor and rice systematically deviate from ECDF.
#' # RMSE is better to sort by than R2
#' 
#' # custom weights
#' cw <- c("gpa"=7, "gev"=3, "wak"=6, "wei"=4, "kap"=3.5, "gum"=3, "ray"=2.1,
#'         "ln3"=2, "pe3"=2.5, "gno"=4, "gam"=5)
#' dlf <- distLfit(annMax, plot=FALSE, weightc=cw)
#' distLgofPlot(dlf, ranks=TRUE)
#' 
#'
#' @param rmse    Numeric: Named vector with goodness of fit values (RMSE).
#'                Can also be a data.frame, in which case the column rmse or RMSE is used.
#' @param order   Logical: should result be ordered by RMSE? DEFAULT: TRUE
#' @param weightc Optional: a named vector with custom weights for each distribution.
#'                Are internally normalized to sum=1 after removing nonfitted dists.
#'                Names match the parameter names from \code{rmse}.
#'                DEFAULT: NA
#'
distLweights <- function(
rmse,
order=TRUE,
weightc=NA
)
{
# get data.frame column:
if(is.data.frame(rmse) | is.matrix(rmse))
  {
  colm <- grep("rmse", colnames(rmse), ignore.case=TRUE)
  if(length(colm)!=1) stop("There is not a single column matching 'rmse' among ", 
                           toString(colnames(rmse)))
  rmse2 <- rmse[,colm]
  names(rmse2) <- rownames(rmse)
  rmse <- rmse2
  }
  
if(is.null(names(rmse))) stop("rmse must have names.")

# the lower RMSE, the better GOF, the more weight
maxrmse <- max(rmse, na.rm=TRUE)
  
# Zero weight to worst fit (needs 2 or more distributions to work):
weight2 <- maxrmse - rmse

# at least a little weight for all distributions:
weight1 <- maxrmse - rmse + min(rmse, na.rm=TRUE)  
# with min or mean added, the worst fit is not completely excluded

# use only best half (needs 3 or more values):
weight3 <-  maxrmse - rmse
weight3[weight3<median(weight3)] <- 0

# custom weight:
if(any(!is.na(weightc)))
  {
  cn <- names(weightc) # custom names
  rn <- names(rmse)
  if(is.null(cn)) stop("weightc must have names.")
  miss <- ! rn %in% cn
  if(any(miss)) warning("names present in rmse, but not in weightc, thus given zero weight: ", 
                        toString(rn[miss]))
  miss <- ! cn %in% rn
  if(any(miss)) 
    {
    warning("names present in weightc, but not in rmse, thus ignored: ", toString(cn[miss]))
    weightc <- weightc[!miss]
    cn <- names(weightc) 
    }
  weightc2 <- rep(0,length(rmse))
  names(weightc2) <- rn
  weightc2[cn] <- weightc
  weightc <- weightc2
} else
weightc <- rep(NA, length(rmse))

# replace NAs with 0
weight1[!is.finite(weight1)] <- 0
weight2[!is.finite(weight2)] <- 0
weight3[!is.finite(weight3)] <- 0
weightc[!is.finite(weightc)] <- 0

# normalize to get sum=1
weight1 <- weight1/sum(weight1) 
weight2 <- weight2/sum(weight2)
weight3 <- weight3/sum(weight3)
weightc <- weightc/sum(weightc)

# output data.frame:
out <- data.frame(weight1, weight2, weight3, weightc)

# order by GOF:
if(order) out <- out[ order(rmse), ] # sorting by R2 does not work, see examples

out
}
