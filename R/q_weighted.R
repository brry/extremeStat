#' Compute weighted averages of quantile estimates
#' 
#' @return data.frame with rows "weighted*" added.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2016
#' @seealso \code{\link{distLquantile}}
#' @keywords distribution
#' @export
#' @examples
#' x <- data.frame(A=1:5, RMSE=runif(5))
#' distLweights(x, onlydn=FALSE)
#' 
#' q_weighted(x,  onlydn=FALSE)
#' q_weighted(x,  distLweights(x, weightc=c("1"=3, "3"=5), order=FALSE, onlydn=FALSE)  )
#' 
#' \dontrun{ # time consuming
#' x <- rexp(190)
#' d <- distLquantile(x)
#' d2 <- q_weighted(d)
#' stopifnot(all(d==d2, na.rm=TRUE))
#' 
#' # fast option for adding custom weighted estimates:
#' cw <- runif(17)
#' names(cw) <- c("exp", "gam", "gev", "glo", "gno", "gpa", "gum", "kap", "lap",
#'                "ln3", "nor", "pe3", "ray", "revgum", "rice", "wak", "wei")
#' dw <- distLweights(d, weightc=cw)
#' qw1 <- q_weighted(d, weightc=cw); qw1
#' qw2 <- q_weighted(d, weights=dw); qw2
#' stopifnot(all(qw1==qw2, na.rm=TRUE))
#' q_weighted(d, weights=dw, onlyc=TRUE)
#' q_weighted(d, weights=data.frame(weightc=cw), onlyc=TRUE)
#' 
#' system.time(pbreplicate(5000, q_weighted(d, weightc=cw)))             # 8.5 secs
#' system.time(pbreplicate(5000, q_weighted(d, weights=dw, onlyc=TRUE))) # 0.8 secs
#' }
#' 
#' @param quant   Data.frame as in \code{\link{distLquantile}} output.
#' @param weights Data.frame as in \code{\link{distLweights}} output.
#' @param onlyc   Logical: only return custom weighted quantile estimates as
#'                a vector? Useful to add those to existing results. See examples.
#'                DEFAULT: FALSE
#' @param \dots   Arguments passed to \code{\link{distLweights}} like
#'                weightc, onlydn=FALSE. order will be ignored, as
#'                q_weighted only adds/changes the rows weighted*.
#' 
q_weighted <- function(
quant,
weights=distLweights(quant, ...),
onlyc=FALSE,
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

if(onlyc) return(Qweighted("weightc"))

quant["weighted1",notRMSE] <- Qweighted("weight1")
quant["weighted2",notRMSE] <- Qweighted("weight2")
quant["weighted3",notRMSE] <- Qweighted("weight3")
quant["weightedc",notRMSE] <- Qweighted("weightc")
quant
}
