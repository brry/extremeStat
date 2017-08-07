#' print dlf objects
#' 
#' print list objects created in this package
#' 
#' The common object to share between functions (see overview in \code{\link{extremeStat}})
#' is a list with the following elements:
#' \tabular{ll}{
#' \code{dat}          \tab numeric vector with (extreme) values,
#'                          with all NAs and values below threshold removed\cr
#' \code{dat_full}     \tab original input data complete with NAs\cr
#' \code{datname}      \tab character string for main, xlab etc \cr
#' \code{parameter}    \tab list (usually of length 17 if \code{speed=TRUE} in
#'                          \code{\link{distLfit}})
#'                          with parameters of each distribution\cr
#' \code{gof}          \tab dataframe with 'Goodness of Fit' measures, sorted by
#'                          RMSE of theoretical and empirical cumulated density\cr
#' \code{distnames}    \tab character vector with selected distribution names\cr
#' \code{distfailed}   \tab Names of nonfitted distributions or ""\cr
#' \code{distcols}     \tab colors for distnames (for plotting). If not given manually,
#'                          determined by \code{berryFunctions::\link{rainbow2}}\cr
#' \code{distselector} \tab character string with function name creating
#'                          the selection\cr
#' \code{truncate, threshold} \tab Truncation percentage and threshold value,
#'                          relevant for \code{\link{distLquantile}}\cr
#' }
#' 
#' optionally, it can also contain:
#' 
#' \tabular{ll}{
#' \code{returnlev, npy }  \tab dataframe with values of distributions for given
#'                          return periods (\code{RPs}), number of observations per year/block.
#'                          These elements are only added in \code{\link{distLextreme}}\cr
#' \code{RPweibull, RPgringorton} \tab Return periods according to plotting positions,
#'                          added in \code{\link{plotLextreme}}\cr
#' \code{quant}        \tab Quantile estimates from \code{\link{distLquantile}}\cr
#' \code{exBootRPs, qexBootSim, exBootCI, exBootCL} \tab objects from \code{\link{distLexBoot}}\cr
#' }
#' 
#' @return none, prints via \code{\link{message}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014, March + July 2015, Dec 2016
#' @seealso \code{\link{extremeStat}}
#' @keywords list methods print
#' @importFrom berryFunctions truncMessage
#' @export
#' @examples
#' 
#' # see
#' ?distLextreme
#' 
#' @param dlf List as explained in section Details
#' @param digits number of digits \code{\link{round}}ed to. DEFAULT: 1
#' 
printL <- function(
dlf,
digits=1
)
{
# input checks:
obj <- deparse(substitute(dlf))
if(!is.list(dlf)) stop("dlf is not a list")
must <- c("dat", "dat_full", "datname", "parameter", "gof",
          "distnames","distcols", "distselector", "distfailed", "truncate", "threshold")
if(any(!must %in% names(dlf))) warning("dlf must include the element(s) ",
                                       toString(must[!must %in% names(dlf)]))
# functions:
vals <- function(x, signi=FALSE, pre="min/median/max: ")
  suppressWarnings({
  nNA <- sum(is.na(x))
  x <- x[!is.na(x)]
  if(!signi)values <-  round(c(min(x), median(x), max(x)), digits)
  if(signi) values <- signif(c(min(x), median(x), max(x)), digits+1)
  paste(pre, paste(values, collapse="/"), " nNA:", nNA)
  })
# message preparation:
n <- length(dlf$dat)
inparnotgof <- ! names(dlf$parameter) %in% rownames(dlf$gof)
ingofnotpar <- ! rownames(dlf$gof) %in% names(dlf$parameter)
PP <- "RPweibull" %in% names(dlf) | "RPgringorton" %in% names(dlf)
RP <- "returnlev" %in% names(dlf)
EB <- "exBootRPs" %in% names(dlf)
QQ <-     "quant" %in% names(dlf)
other <- ! names(dlf) %in% c(must, "returnlev", "RPweibull", "RPgringorton",
                             "exBootRPs", "exBootSim", "exBootCI", "exBootCL",
                             "quant", "npy")
if(RP) rlev <- substr(colnames(dlf$returnlev), start=4, stop=8)

# message output:
message("----------\nDataset '",dlf$datname,"' with ",n," values. ",vals(dlf$dat),
if( ! is.vector(dlf$dat)) paste0("\n--> dat is not a vector, but ",class(dlf$dat),"!"),
"\ntruncate: ", dlf$truncate, " threshold: ",round(dlf$threshold,digits+1),
    ". dat_full with ", length(dlf$dat_full), " values: ", vals(dlf$dat_full, pre=""),
"\ndlf with ", nrow(dlf$gof), " distributions. In descending order of fit quality:\n",
   toString(rownames(dlf$gof)[order(dlf$gof$RMSE)]),
if(any(inparnotgof)) paste0("\n--> dists in parameter but not in gof: ",
                      toString(names(dlf$parameter)[inparnotgof]) ),
if(any(ingofnotpar)) paste0("\n--> dists in gof but not in parameter: ",
                       toString(rownames(dlf$gof)[ingofnotpar])),
"\nRMSE ", vals(dlf$gof$RMSE, TRUE),
"\n", length(dlf$distnames), " distnames + ", length(dlf$distcols),
      " distcols from distselector ", dlf$distselector,
if(dlf$distfailed[1]!="") paste0("\n--> fitting failed for ",length(dlf$distfailed),
                                " distributions: ", toString(dlf$distfailed)),
if(QQ) paste("\nquant:",nrow(dlf$quant),"rows,",ncol(dlf$quant),"columns,",
            prod(dim(dlf$quant)),"values, of which",sum(is.na(dlf$quant)),"NA."),
if(PP | RP) "\n-- distLextreme elements:",
if(PP) paste0("\n",length(dlf$RPweibull)," Plotting positions ",
              vals(c(dlf$RPweibull, dlf$RPgringorton))),
if(RP) paste0("\n", ncol(dlf$returnlev), " Return Periods (", vals(as.numeric(rlev)), "): ",
              berryFunctions::truncMessage(rlev, ntrunc=10, prefix="", midfix=""),
              "\nReturn Levels ", vals(head(dlf$returnlev,-3)), "\nnpy: ", dlf$npy),
if(EB) paste0("\n-- distLexBoot elements (",ncol(dlf$exBootSim[[1]])," simulations @ ",
              dlf$exBootCL*100, "% conf.lev):\n",
              length(dlf$exBootSim)," distributions: ", toString(names(dlf$exBootSim)),
              "\n", length(dlf$exBootRPs), " Return periods ", vals(dlf$exBootRPs),
              "\nReturn levels ", vals(unlist(dlf$exBootSim))),
if(any(other)) paste0("\n--> Other elements in the object '", obj, "': ",
                      toString(names(dlf)[other]))
)
}
