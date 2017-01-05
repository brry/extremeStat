#' Compute weighted averages of quantile estimates
#'
#' @return data.frame with rows "weighted*" added.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLquantile}}
#' @keywords distribution
#' @export
#' @examples
#' x <- data.frame(A=1:5, RMSE=runif(5))
#' distLweights(x)
#' 
#' q_weighted(x)
#' q_weighted(x,  distLweights(x, weightc=c("1"=3, "3"=5), order=FALSE)  )
#' 
#' x <- rexp(190)
#' d <- distLquantile(x)
#' d2 <- q_weighted(d)
#' stopifnot(all(d==d2, na.rm=TRUE))
#'
#' @param quant Data.frame as in \code{\link{distLquantile}} output.
#' @param weights Data.frame as in \code{\link{distLweights}} output.
#' @param \dots   Arguments passed to \code{\link{distLweights}} 
#'                like weightc, onlydn=FALSE. order will be ignored, as 
#'                q_weighted only adds/changes the rows weighted*.
#' 
q_weighted <- function(
quant,
weights=distLweights(quant, ...),
...
)
{
weights <- weights # evaluate promise
notRMSE <- which(colnames(quant)!="RMSE")
Qweighted <- function(w)
  {
  sapply(notRMSE, function(col_n)
    {
    vals <- quant[rownames(weights),col_n]
    wx <- weights[,w]
    finite <- is.finite(vals)  & is.finite(wx)
    if(!any(finite)) return(NA)
    wx[!finite] <- 0
    vals[!finite] <- 0
    wx <- wx/sum(wx) # rescale to 1
    sum(vals*wx)
    })
  }
quant["weighted1",notRMSE] <- Qweighted("weight1")
quant["weighted2",notRMSE] <- Qweighted("weight2")
quant["weighted3",notRMSE] <- Qweighted("weight3")
quant["weightedc",notRMSE] <- Qweighted("weightc")
quant
}
