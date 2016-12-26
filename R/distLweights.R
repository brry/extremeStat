#' Compute distribution weights from GOF
#'
#' Determine distribution function weights from RMSE for weighted averages.
#' The weights are inverse to RMSE: weight1 for all dists, 
#' weight2 places zero weight on the worst fitting function, 
#' weight3 on the worst half of functions.
#'
#' @return data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLgof}}, \code{\link{distLquantile}}
#' @keywords distribution
#' @export
#' @examples
#' # weights from RMSE vector:
#' rmse <- c(gum=0.20, wak=0.17, gam=0.21, gev=0.15)
#' distLweights(rmse)
#' distLweights(rmse, order=FALSE)
#' 
#' # weights from RMSE in data.frame:
#' df <- data.frame("99.9%"=2:5, rmse=sample(3:6))
#' rownames(df) <- letters[1:4]
#' df ;  distLweights(df)
#' 
#' # custom weights:
#' set.seed(42); x <- data.frame(A=1:5, RMSE=runif(5)) ; x
#' distLweights(x)
#' distLweights(x, weightc=c("1"=3, "3"=5)) 
#' distLweights(x, weightc=c("1"=3, "3"=5), order=FALSE) 
#' 
#' # real life example:
#' data(annMax)
#' cw <- c("gpa"=7, "gev"=3, "wak"=6, "wei"=4, "kap"=3.5, "gum"=3, "ray"=2.1,
#'         "ln3"=2, "pe3"=2.5, "gno"=4, "gam"=5)
#' dlf <- distLfit(annMax, plot=FALSE, weightc=cw, quiet=TRUE, order=FALSE)
#' plotLgof(dlf)
#' 
#' # GOF judgement by RMSE, not R2 --------
#' # Both RMSE and R2 are computed with ECDF and TCDF
#' # R2 may be very good (see below), but fit needs to be close to 1:1 line, 
#' # which is better measured by RMSE
#' set.seed(42); x <- runif(30)
#' y1 <-     x+rnorm(30,sd=0.09)
#' y2 <- 1.5*x+rnorm(30,sd=0.01)-0.3
#' plot(x,x, asp=1, las=1)
#' berryFunctions::linReg(x, y2, add=TRUE, digits=4, pos1="topleft")
#' points(x,y2, col="red", pch=3)
#' points(x,y1, col="blue")
#' berryFunctions::linReg(x, y1, add=TRUE, digits=4, col="blue", pos1="left")
#' abline(a=0, b=1, lwd=3, lty=3)
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
  if(length(colm)<1) stop("There is no column matching 'rmse' among ", 
                           toString(colnames(rmse)))
  if(length(colm)>1) stop("There are several columns matching 'rmse' among ", 
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
if(sum(weight2,na.rm=TRUE)==0) weight2[weight2==0] <- 1

# at least a little weight for all distributions:
weight1 <- maxrmse - rmse + min(rmse, na.rm=TRUE)  
# with min or mean added, the worst fit is not completely excluded
if(sum(weight1,na.rm=TRUE)==0) weight1[weight1==0] <- 1

# use only best half (needs 3 or more values):
weight3 <-  weight2
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
mysum <- function(x) {y <- sum(x); if(y==0) 1 else y}
weight1 <- weight1/mysum(weight1) 
weight2 <- weight2/mysum(weight2)
weight3 <- weight3/mysum(weight3)
weightc <- weightc/mysum(weightc)

# output data.frame:
out <- data.frame(rmse, weight1, weight2, weight3, weightc)

# order by GOF:
if(order) out <- out[ order(rmse), ] # sorting by R2 does not work, see examples

out
}
