#' Deprecated Functions
#' 
#' With version 1.3.0 (Jan 2017), a number of functions have been renamed and restructured.
#' The computing functions don't plot anymore. 
#' This is changed to reduce number of arguments in functions 
#' and increase reproducibility by explicitly creating dlf objects.
#' Renaming has been done to create a consistent and short naming pattern.
#' See the new structure in \code{\link{extremeStat}}.
#'
#' @name extremeStat-deprecated
#' @aliases distLextremePlot distLgofPlot distLplot distLprint
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2017
#' @seealso \code{\link{extremeStat}}
#' @param \dots Ignored arguments
#' @export
distLextremePlot <- function(...) deprmes("plotLextreme")
#' @export
#' @rdname extremeStat-deprecated
distLgof <- function(...) deprmes("distLfit")
#' @export
#' @rdname extremeStat-deprecated
distLgofPlot <- function(...) deprmes("plotLweights")
#' @export
#' @rdname extremeStat-deprecated
distLplot <- function(...) deprmes("plotLfit")
#' @export
#' @rdname extremeStat-deprecated
distLprint <- function(...) deprmes("printL")

deprmes <- function(new) warning("'", as.character(sys.call(sys.parent()))[1L],
                              "' is deprecated. Instead, use  ", new, "\n",
                              "help('extremeStat-deprecated')", call.=FALSE)


if(FALSE){
".onAttach" <- function(lib, pkg) 
{
  d <- utils::packageDescription(pkg)
  packageStartupMessage("# Loaded extremeStat ", d$Version, " (",d$Date,"). Package restructured since 0.6.0 (2016-12-13).\n",
                        "# Computing functions don't plot anymore and some are renamed. See help('extremeStat-deprecated')")
}
}

