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
#' @param quant Data.frame as in \code{\link{distLquantile}} output.
#' @param weights Data.frame as in \code{\link{distLweights}} output.
#'
q_weighted <- function(
quant,
weights=distLweights(quant)
)
{
weights <- weights # evaluate promise
Qweighted <- function(w)
  {
  sapply(1:ncol(quant), function(col_n)
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
quant["weighted1",] <- Qweighted("weight1")
quant["weighted2",] <- Qweighted("weight2")
quant["weighted3",] <- Qweighted("weight3")
quant["weightedc",] <- Qweighted("weightc")
quant
}
