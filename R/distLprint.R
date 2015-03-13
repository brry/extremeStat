# print dlf objects
# Berry Boessenkool, Sept 2014

distLprint <- function(
dlf, # List as explained in \code{\link{extremeStat}}
digits=1 # number of digits \code{\link{round}}ed to.
)
{
# input checks:
obj <- deparse(substitute(dlf))
if(!is.list(dlf)) stop("dlf is not a list")
must <- c("dat", "datname",  "gofProp", "parameter", "gof")
for(i in 1:5) if(!must[i] %in% names(dlf)) stop("dlf must include the element '", must[i], "'.")
# functions:
pastec <- function(...) paste(..., collapse=", ")
vals <- function(x, signi=FALSE) 
  {
  nNA <- sum(is.na(x))
  x <- x[!is.na(x)]
  if(!signi)values <-  round(c(min(x), median(x), max(x)), digits)
  if(signi) values <- signif(c(min(x), median(x), max(x)), digits+1)
  paste(pastec(values), " nNA:", nNA)
  } 
# message preparation:
n <- length(dlf$dat)
inparnotgof <- ! names(dlf$parameter) %in% rownames(dlf$gof)
ingofnotpar <- ! rownames(dlf$gof) %in% names(dlf$parameter)
PP <- "RPweibull" %in% names(dlf) | "RPgringorton" %in% names(dlf)
RP <- "returnlev" %in% names(dlf)
other <- ! names(dlf) %in% c(must, "RPweibull", "RPgringorton" , "returnlev")
# Further elements:
if(PP) PPs <- vals(c(dlf$RPweibull, dlf$RPgringorton))
if(RP) RPs <- substr(colnames(dlf$returnlev), start=4, stop=8)

# message output:
message("Dataset '", dlf$datname, "' with ", n, " values. min/median/max: ", vals(dlf$dat),
if( ! is.vector(dlf$dat)) "\n--> dat is not a vector!",
"\ndlf with ", nrow(dlf$gof), " distributions. In descending order of fit quality:\n", 
   pastec(rownames(dlf$gof)),
if(any(inparnotgof)) "\n--> dists in parameter but not in gof: ",
if(any(inparnotgof)) paste(names(dlf$parameter)[inparnotgof], collapse=", "),
if(any(ingofnotpar)) "\n--> dists in gof but not in parameter: ",
if(any(ingofnotpar)) paste(rownames(dlf$gof)[ingofnotpar], collapse=", "),
"\ngofProp: ", dlf$gofProp, ", RMSE min/median/max: ", vals(dlf$gof$RMSE, TRUE), 
if(PP | RP) "\ndistLextreme elements:",
if(PP) "\nPlotting positions min/median/max: ",
if(PP) PPs,
if(RP) "\nReturn Periods: ",
if(RP) pastec(RPs),
if(RP) "\nReturn Levels min/median/max: ",
if(RP) vals(dlf$returnlev),
if(any(other)) paste("\n--> Other elements in the object '", obj, "': "),
if(any(other)) paste(names(dlf)[other], collapse=", ")
)
}
