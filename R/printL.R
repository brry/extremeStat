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
#' \code{returnlev, }  \tab dataframe with values of distributions for given 
#'                          return periods (\code{RPs}). This element is only added
#'                          in \code{\link{distLextreme}}\cr
#' \code{RPweibull, RPgringorton} \tab Return periods according to plotting positions, 
#'                          added in \code{\link{plotLextreme}}\cr
#' \code{quant}        \tab Quantile estimates from \code{\link{distLquantile}}\cr
#' \code{exBootRPs, qexBootSim, exBootCI} \tab objects from \code{\link{distLexBoot}}\cr
#' }
#' 
#' @return none, prints via \code{\link{message}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014, March + July 2015, Dec 2016
#' @seealso \code{\link{extremeStat}}
#' @keywords list methods print
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
          "distnames","distcols", "distselector", "truncate", "threshold")
if(any(!must %in% names(dlf))) warning("dlf must include the element(s) ", 
                                       toString(must[!must %in% names(dlf)]))
# functions:
vals <- function(x, signi=FALSE)
  suppressWarnings({
  nNA <- sum(is.na(x))
  x <- x[!is.na(x)]
  if(!signi)values <-  round(c(min(x), median(x), max(x)), digits)
  if(signi) values <- signif(c(min(x), median(x), max(x)), digits+1)
  paste(paste(values, collapse="/"), " nNA:", nNA)
  })
# message preparation:
n <- length(dlf$dat)
inparnotgof <- ! names(dlf$parameter) %in% rownames(dlf$gof)
ingofnotpar <- ! rownames(dlf$gof) %in% names(dlf$parameter)
PP <- "RPweibull" %in% names(dlf) | "RPgringorton" %in% names(dlf)
RP <- "returnlev" %in% names(dlf)
qq <- "quant" %in% names(dlf)
other <- ! names(dlf) %in% c(must, "returnlev", "RPweibull", "RPgringorton", 
                             "quant", "exBootRPs", "qexBootSim", "exBootCI")


# message output:
message("----------\nDataset '", dlf$datname, "' with ", 
        n, " values. min/median/max: ", vals(dlf$dat),
if( ! is.vector(dlf$dat)) "\n--> dat is not a vector!",
"\ntruncate: ", dlf$truncate, " threshold: ",round(dlf$threshold,digits+1),
    ". dat_full with ", length(dlf$dat_full), " values: ", vals(dlf$dat_full),
"\ndlf with ", nrow(dlf$gof), " distributions. In descending order of fit quality:\n", 
   toString(rownames(dlf$gof)[order(dlf$gof$RMSE)]),
if(any(inparnotgof)) paste0("\n--> dists in parameter but not in gof: ",
                      toString(names(dlf$parameter)[inparnotgof]) ),
if(any(ingofnotpar)) paste0("\n--> dists in gof but not in parameter: ",
                       toString(rownames(dlf$gof)[ingofnotpar])),
"\nRMSE min/median/max: ", vals(dlf$gof$RMSE, TRUE),
"\n", length(dlf$distnames), " distnames + ", length(dlf$distcols), 
      " distcols from distselector ", dlf$distselector,
if(qq) paste("\nquant:",nrow(dlf$quant),"rows,",ncol(dlf$quant),"columns,",
            prod(dim(dlf$quant)),"values, of which",sum(is.na(dlf$quant)),"NA."),
if(PP | RP) "\ndistLextreme elements:",
if(PP) paste0("\nPlotting positions min/median/max: ", vals(c(dlf$RPweibull, dlf$RPgringorton))),
if(RP) paste0("\nReturn Periods: ", toString(substr(colnames(dlf$returnlev), start=4, stop=8)),
              "\nReturn Levels min/median/max: ", vals(dlf$returnlev)),
if(any(other)) paste0("\n--> Other elements in the object '", obj, "': ",
                      toString(names(dlf)[other]))
)
}
