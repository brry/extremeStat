#' Compute distribution weights from GOF
#'
#' Determine distribution function weights from RMSE for weighted averages.
#'
#' @return data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{help}}, \code{\link{help}}
#' @keywords distribution
#' @export
#' @examples
#' 
#' distLweights(c(gum=0.20, wak=0.17, gam=0.21))
#' distLweights(c(gum=0.20, wak=0.17, gam=0.21), order=FALSE)
#' distLweights(c(gum=0.20))
#' 
#' # RMSE vs R2 for GOF judgement --------
#' library(lmomco)
#' data(annMax)
#' 
#' compranks <- function(d)
#' {
#' gofProp <- 0.5
#' x <- sort(annMax, decreasing=TRUE)[  1:(gofProp*length(annMax))  ]
#' tcdfs <- plmomco(x,dlf50$parameter[[d]])
#' ecdfs <- ecdf(annMax)(x) # Empirical CDF
#' # Root Mean Square Error, R squared:
#' berryFunctions::linReg(tcdfs, ecdfs, lwd=1, pch=16, main=d, digits=5, xlim=c(0.5, 1),
#'        ylim=c(0.5, 1), pos1="topleft")
#' abline(a=0,b=1, lty=3)
#' c(berryFunctions::rmse(tcdfs, ecdfs), berryFunctions::rsquare(tcdfs, ecdfs))
#' }
#' 
#' dlf50 <- distLfit(annMax, gofProp=0.5)
#' dn <- rownames(dlf50$gof)
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
#' @param rmse  Numeric: (named) vector with goodness of fit values (RMSE)
#' @param order Logical: should result be ordered by RMSE? DEFAULT: TRUE
#' @param weightc Numeric: Named vector with custom weights. DEFAULT: NA
#'
distLweights <- function(
rmse,
order=TRUE,
weightc=NA
)
{

# the lower RMSE, the better GOF, the more weight
  
if(is.null(names(rmse))) stop("rmse must have names.")
  
# Zero weight to worst fit (needs 2 or more distributions to work):
weight2 <- rmse
weight2 <- max(weight2) - weight2

# at least a little weight for all distributions:
weight1 <- rmse 
weight1 <- max(weight1) - weight1 + min(weight1)  
# with min or mean added, the worst fit is not completely excluded

# use only best half (needs 4 or more values):
weight3 <- rmse
weight3 <- max(weight3) - weight3
n <- length(rmse)
weight3[(n/2):n] <- 0

# custom weight:
if(any(!is.na(weightc)))
  {
  cn <- names(weightc) # custom names
  rn <- names(rmse)
  if(is.null(cn)) stop("weightc must have names.")
  miss <- ! rn %in% cn
  if(any(miss)) warning("names present in rmse, but not in weightc, thus given weight 0: ", 
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

# normalize to get sum=1
weight1 <- weight1/sum(weight1) 
weight2 <- weight2/sum(weight2)
weight3 <- weight3/sum(weight3)
weightc <- weightc/sum(weightc)

# output data.frame:
out <- data.frame(weight1, weight2, weight3, weightc)
#rownames(out) <- names(rmse)

# order by GOF:
if(order) out <- out[ order(rmse), ] # sorting by R2 does not work, see examples

out
}
