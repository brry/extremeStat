#' print dlf objects
#' 
#' print list objects created in this package
#' 
#' @return none, prints via \code{\link{message}}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014, March + July 2015
#' @seealso \code{\link{extremeStat}}
#' @keywords list methods print
#' @export
#' @examples
#' 
#' # see 
#' ?distLextreme
#' 
#' @param dlf List as explained in \code{\link{extremeStat}}
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
must <- c("dat", "dat_full", "datname", "parameter", "gof", "truncate", "threshold")
for(i in 1:length(must)) if(!must[i] %in% names(dlf)) warning("dlf must include the element '", must[i], "'.")
# functions:
vals <- function(x, signi=FALSE)
  {
  nNA <- sum(is.na(x))
  x <- x[!is.na(x)]
  if(!signi)values <-  round(c(min(x), median(x), max(x)), digits)
  if(signi) values <- signif(c(min(x), median(x), max(x)), digits+1)
  paste(paste(values, collapse="/"), " nNA:", nNA)
  } 
# message preparation:
n <- length(dlf$dat)
inparnotgof <- ! names(dlf$parameter) %in% rownames(dlf$gof)
ingofnotpar <- ! rownames(dlf$gof) %in% names(dlf$parameter)
PP <- "RPweibull" %in% names(dlf) | "RPgringorton" %in% names(dlf)
RP <- "returnlev" %in% names(dlf)
CD <- "coldist" %in% names(dlf)
qq <- "quant" %in% names(dlf)
other <- ! names(dlf) %in% c(must, "RPweibull", "RPgringorton" , "returnlev", "coldist", "quant")
# Further elements:
if(PP) PPs <- vals(c(dlf$RPweibull, dlf$RPgringorton))
if(RP) RPs <- substr(colnames(dlf$returnlev), start=4, stop=8)

# message output:
message("----------\nDataset '", dlf$datname, "' with ", n, " values. min/median/max: ", vals(dlf$dat),
if( ! is.vector(dlf$dat)) "\n--> dat is not a vector!",
"\ntruncate: ", dlf$truncate, " threshold: ",round(dlf$threshold,digits+1),
    ". dat_full with ", length(dlf$dat_full), " values: ", vals(dlf$dat_full),
"\ndlf with ", nrow(dlf$gof), " distributions. In descending order of fit quality:\n", 
   toString(rownames(dlf$gof)),
if(any(inparnotgof)) "\n--> dists in parameter but not in gof: ",
if(any(inparnotgof)) toString(names(dlf$parameter)[inparnotgof]),
if(any(ingofnotpar)) "\n--> dists in gof but not in parameter: ",
if(any(ingofnotpar)) toString(rownames(dlf$gof)[ingofnotpar]),
"\nRMSE min/median/max: ", vals(dlf$gof$RMSE, TRUE),
if(qq) paste("\nquant:",nrow(dlf$quant),"rows,",ncol(dlf$quant),"columns,",
            prod(dim(dlf$quant)),"values, of which",sum(is.na(dlf$quant)),"NA."),
if(CD) paste("\n", length(dlf$coldist),"distribution colors:", toString(dlf$coldist)),
if(PP | RP) "\ndistLextreme elements:",
if(PP) "\nPlotting positions min/median/max: ",
if(PP) PPs,
if(RP) "\nReturn Periods: ",
if(RP) toString(RPs),
if(RP) "\nReturn Levels min/median/max: ",
if(RP) vals(dlf$returnlev),
if(any(other)) paste0("\n--> Other elements in the object '", obj, "': "),
if(any(other)) toString(names(dlf)[other])
)
}
