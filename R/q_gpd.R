#' GPD quantile of sample
#' 
#' Compute quantile of General Pareto Distribution fitted to sample by peak over treshold (POT) method
#' using treshold from truncation proportion,
#' comparing several R packages doing this
#'
#' @details Depending on the value of "package", this fits the GPD using \cr
#' \code{evir::\link[evir]{gpd}}\cr
#' \code{evd::\link[evd]{fpot}}\cr
#' \code{extRemes::\link[extRemes]{fevd}}\cr
#' \code{fExtremes::\link[fExtremes]{gpdFit}}\cr
#' \code{ismev::\link[ismev]{gpd.fit}}\cr
#' \code{Renext::\link[Renext]{Renouv}} or \code{Renext::\link[Renext]{fGPD}}\cr\cr
#'
#' The \bold{\code{method}} defaults (and other possibilities) are \cr
#' evir: "pwm" (probability-weighted moments), or "ml" (maximum likelihood) \cr
#' evd: none, only Maximum-likelihood fitting implemented \cr
#' extRemes: "MLE", or "GMLE", "Bayesian", "Lmoments" \cr
#' fExtremes: "pwm", or "mle"\cr
#' ismev: none, only Maximum-likelihood fitting implemented \cr
#' Renext: "r" for \code{\link[Renext]{Renouv}}, or 'f' (no truncation, all negative values ignored!) for \code{\link[Renext]{fGPD}}\cr\cr
#'
#' The Quantiles are always given with \code{probs} in regard to the full (uncensored) sample.
#' If e.g. truncate is 0.90, the distribution function is fitted to the top 10\% of the sample.
#' The 95th percentile of the full sample is equivalent to the 50\% quantile of the subsample acutally used for fitting.
#' For computation, the probabilities are internally updated with \code{p2=(p-t)/(1-t)}
#' but labelled with the original \code{p}.
#' If you truncate 90\% of the sample, you cannot compute the 70th percentile anymore,
#' thus \code{undertruncNA} should be left to TRUE. \cr
#' If not exported by the packages, the quantile functions are extracted from their current (Feb 2016) source code.
#' 
#' @return Named vector of quantile estimates for each value of \code{probs},\cr
#'    or if(returnlist): list with element \code{q_gpd_quant} and info-elements added.
#'    q_gpd_n_geq is number of values greater than or equal to q_gpd_threshold. gt is only greater than.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @references \url{http://stackoverflow.com/questions/27524131/calculation-of-return-levels-based-on-a-gpd-in-different-r-packages}\cr
#'             \url{http://stats.stackexchange.com/questions/129438/different-quantiles-of-a-fitted-gpd-in-different-r-packages}
#' @seealso \code{\link{distLquantile}} which compares results for all packages\cr
#'          Other related packages (not implemented):\cr
#'                \url{https://cran.r-project.org/package=gPdtest}\cr
#'                \url{https://cran.r-project.org/package=actuar}\cr
#'                \url{https://cran.r-project.org/package=fitdistrplus}\cr
#'                \url{https://cran.r-project.org/package=lmom}\cr
#' @keywords distribution robust univar
#' @export
#' @importFrom evir quant gpd
#' @importFrom evd fpot qgpd
#' @importFrom extRemes fevd qevd
#' @importFrom fExtremes gpdFit
#' @importFrom ismev gpd.fit
#' @importFrom Renext fGPD Renouv qGPD
#' @importFrom berryFunctions quantileMean
#'
#' @examples
#' data(annMax)
#' q_gpd(annMax)
#' q_gpd(annMax, truncate=0.6)
#' q_gpd(annMax, truncate=0.85)
#' q_gpd(annMax, truncate=0.91) 
#'
#' q_gpd(annMax, package="evir")
#' q_gpd(annMax, package="evir", method="ml")
#' q_gpd(annMax, package="evd")
#' q_gpd(annMax, package="extRemes")
#' q_gpd(annMax, package="extRemes", method="GMLE")
#' #q_gpd(annMax, package="extRemes", method="Bayesian") # computes a while
#' q_gpd(annMax, package="extRemes", method="Lmoments")
#' q_gpd(annMax, package="extRemes", method="nonsense") # NAs
#' q_gpd(annMax, package="fExtremes")                   # log warnings
#' q_gpd(annMax, package="fExtremes", efquiet=TRUE)    # silenced warnings
#' q_gpd(annMax, package="fExtremes", method= "mle")
#' q_gpd(annMax, package="ismev")
#' q_gpd(annMax, package="Renext")
#' q_gpd(annMax, package="Renext", method="f")
#' dummy <- try(q_gpd(annMax, package="nonsense"), silent=TRUE) # error
#' stopifnot(class(dummy)=="try-error")
#'
#' q_gpd(annMax, truncate=0.85, package="evd")          # Note about quantiles
#' q_gpd(annMax, truncate=0.85, package="evir")
#' q_gpd(annMax, truncate=0.85, package="evir", quiet=TRUE) # No note
#' q_gpd(annMax, truncate=0.85, package="evir", undertruncNA=FALSE)
#'
#' q_gpd(annMax, truncate=0.85, package="evir", returnlist=TRUE)
#' str(  q_gpd(annMax, truncate=0.85, probs=0.6, package="evir", returnlist=TRUE) )# NAs
#' str(  q_gpd(annMax, package="evir",      returnlist=TRUE)   )
#' str(  q_gpd(annMax, package="evd",       returnlist=TRUE)   )
#' str(  q_gpd(annMax, package="extRemes",  returnlist=TRUE)   )
#' str(  q_gpd(annMax, package="fExtremes", returnlist=TRUE)   )
#' str(  q_gpd(annMax, package="ismev",     returnlist=TRUE)   )
#' str(  q_gpd(annMax, package="Renext",    returnlist=TRUE)   )
#'
#' q_gpd(annMax, package="evir", truncate=0.9, method="ml") # NAs (MLE fails often)
#'
#'
#' \dontrun{
#' ## Not run in checks because simulation takes too long
#' # truncation effect
#' mytruncs <- seq(0, 0.9, len=150)
#' oo <- options(show.error.messages=FALSE, warn=-1)
#' myquants <- sapply(mytruncs, function(t) q_gpd(annMax, truncate=t, quiet=TRUE))
#' options(oo)
#' plot(1, type="n", ylim=range(myquants, na.rm=TRUE), xlim=c(0,0.9), las=1,
#'      xlab="truncated proportion", ylab="estimated quantiles")
#' abline(h=quantileMean(annMax, probs=c(0.8,0.9,0.99)))
#' for(i in 1:3) lines(mytruncs, myquants[i,], col=i)
#' text(0.3, c(87,97,116), rownames(myquants), col=1:3)
#'
#'
#' # Underestimation in small samples
#' # create known population:
#' dat <- extRemes::revd(1e5, scale=50, shape=-0.02, threshold=30, type="GP")
#' op <- par(mfrow=c(1,2), mar=c(2,2,1,1))
#' hist(dat, breaks=50, col="tan")
#' berryFunctions::logHist(dat, breaks=50, col="tan")
#' par(op)
#'
#' # function to estimate empirical and GPD quantiles from subsamples
#' samsizeeffect <- function(n, nrep=30, probs=0.999, trunc=0.5, Q=c(0.4,0.5,0.6))
#' {
#' res <- replicate(nrep, {
#' subsample <- sample(dat, n)
#' qGPD <- q_gpd(subsample, probs=probs, truncate=trunc)
#' qEMP <- berryFunctions::quantileMean(subsample, probs=probs, truncate=trunc)
#' c(qGPD=qGPD, qEMP=qEMP)})
#' apply(res, MARGIN=1, berryFunctions::quantileMean, probs=Q)
#' }
#'
#' # Run and plot simulations
#' samplesize <- c(seq(20, 150, 10), seq(200,800, 100))
#' results <- pbapply::pblapply(samplesize, samsizeeffect)
#' res <- function(row, col) sapply(results, function(x) x[row,col])
#' berryFunctions::ciBand(yu=res(3,1),yl=res(1,1),ym=res(2,1),x=samplesize,
#'   main="99.9% Quantile underestimation", xlab="subsample size", ylim=c(200,400), colm=4)
#' berryFunctions::ciBand(yu=res(3,2),yl=res(1,2),ym=res(2,2),x=samplesize, add=TRUE)
#' abline(h=berryFunctions::quantileMean(dat, probs=0.999))
#' text(300, 360, "empirical quantile of full sample")
#' text(300, 340, "GPD parametric estimate", col=4)
#' text(300, 300, "empirical quantile estimate", col="green3")
#'
#' } # end of dontrun
#'
#' @param x Vector with numeric values. NAs are silently ignored.
#' @param probs Probabilities of truncated (Peak over treshold) quantile. DEFAULT: c(0.8,0.9,0.99)
#' @param truncate Truncation percentage (proportion of sample discarded). DEFAULT: 0
#' @param threshold POT cutoff value. If you want correct percentiles, set this only via truncate, see Details. DEFAULT: \code{\link[berryFunctions]{quantileMean}(x, truncate)}
#' @param package Character string naming package to be used. One of c("evir","evd","extRemes","fExtremes","ismev"). DEFAULT: "extRemes"
#' @param method \code{method} passed to the fitting function, if applicable. Defaults are internally specified (See Details), depending on \code{package}, if left to the DEFAULT: NULL.
#' @param returnlist Return result from the fitting funtion with the quantiles added to the list as element \code{quant} and some information in elements starting with \code{q_gpd_}. DEFAULT: FALSE
#' @param undertruncNA Return NAs for probs below truncate? Highly recommended to leave this at the DEFAULT: TRUE
#' @param quiet Should messages from this function be suppressed? DEFAULT: FALSE
#' @param ttquiet Should truncation!=threshold messages from this function be suppressed? DEFAULT: quiet
#' @param efquiet Should warnings in function calls to the external packages be suppressed via \code{\link{options}(warn=-1)}? The usual type of warning is: NAs produced in log(...). DEFAULT: quiet
#' @param \dots Further arguments passed to the fitting funtion listed in section Details.
#'
q_gpd <- function(
x,
probs=c(0.8,0.9,0.99),
truncate=0,
threshold=berryFunctions::quantileMean(x, truncate),
package="extRemes",
method=NULL,
returnlist=FALSE,
undertruncNA=TRUE,
quiet=FALSE,
ttquiet=quiet,
efquiet=quiet,
...)
{
# Input control: ---------------------------------------------------------------
if(length(package)!=1) stop("package must have length 1, not ", length(package))
pospack <- c("evir","evd","extRemes","fExtremes","ismev","Renext")
if(!package %in% pospack) stop("package ('",
     package, "') must be one of:\n  '", paste(pospack, collapse="', '"), "'.")
x <- x[!is.na(x)]
if(length(truncate)>1)
  {
  truncate <- truncate[1] #
  if(!quiet) on.exit(message("Note in q_gpd: only first value of 'truncate' is used."), add=TRUE)
  }
if(truncate>1 | truncate<0) stop("truncate (proportion discarded) must be 0<t<1, not ", truncate)
#
# truncation probs update: -----------------------------------------------------
probs2 <- probs
if(truncate!=0)
  {
  probs2 <- (probs-truncate)/(1-truncate) # correct probabilities for smaller sample proportion
  probs2[probs < truncate] <- 0   # avoid negative values
  }
#
# threshold warning:
if(!missing(threshold))
  {
  normalthr <- berryFunctions::quantileMean(x, truncate)
  if(threshold != normalthr)
    {
    probs2 <- probs
    if(!ttquiet) on.exit(message("Note in q_gpd: threshold (",threshold,
    ") is not equal to threshold computed from truncate (",normalthr,
    ").\n  Probabilities are not corrected for truncation!"), add=TRUE)
    }
  }
#
# prepare failure output + regular output for list:
failout <- rep(NA, length(probs))
names(failout) <- paste0(probs*100,"%")
outlist <- list(q_gpd_created=Sys.time(), q_gpd_creator="",
   q_gpd_truncate=truncate, q_gpd_threshold=threshold, q_gpd_n_full=length(x),
   q_gpd_n_geq=sum(x>=threshold), q_gpd_n_gt=sum(x>threshold),
   q_gpd_probs=probs, q_gpd_probs2=probs2,
   q_gpd_package=package, q_gpd_method=method, q_gpd_quant=failout)
if(returnlist) failout <- c(list(z="not fitted"), outlist)
# check probs for compliance with truncate:
if(all(probs < truncate) & undertruncNA)
   {
   if(!quiet) on.exit(message("Note in q_gpd: With undertruncNA=TRUE, 'probs' (",
              toString(probs),
              ")\n  must contain values that are larger than 'truncate' (",
              truncate, "). Returning NAs."), add=TRUE)
   return(failout)
   }
#
#
# fitting: ---------------------------------------------------------------------
if(efquiet) {oo <- options(warn=-1); on.exit(options(oo), add=TRUE)}
# function to perform in case of failure. Yields useful error message and returns NAs:
failfun <- function(z, fitfun) {
  if(!quiet) on.exit(message(
     "Note in q_gpd: ",fitfun," failed. Returning NAs. \n  Reason: ", z), add=TRUE)
  return(failout)
  }
# actual fitting:
if(package=="evir") ##################
{
  outlist$q_gpd_creator <- "evir::gpd"
  if(is.null(method)) method <- "pwm"
  pos <- sum(x > threshold)
  z <- try(evir::gpd(x, nextremes=pos, method=method, ...), silent=TRUE)
 if(inherits(z, "try-error")) return(failfun(z, "evir::gpd"))
} else
if(package=="evd") ##################
{
  outlist$q_gpd_creator <- "evd::fpot"
  z <- try(evd::fpot(x, threshold=threshold, model="gpd", std.err=FALSE, ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "evd::fpot"))
} else
if(package=="extRemes") ##################
{
  outlist$q_gpd_creator <- "extRemes::fevd"
  if(is.null(method)) method <- "MLE"
  z <- try(extRemes::fevd(x, method=method, type="GP", threshold=threshold, ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "extRemes::fevd"))
} else
if(package=="fExtremes") ##################
{
  outlist$q_gpd_creator <- "fExtremes::gpdFit"
  if(is.null(method)) method <- "pwm"
  z <- try(z <- fExtremes::gpdFit(x, type=method, u=threshold, ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "fExtremes::gpdFit"))
} else
if(package=="ismev") ##################
{
  outlist$q_gpd_creator <- "ismev::gpd.fit"
  z <- try(ismev::gpd.fit(x, threshold=threshold, show=FALSE, ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "ismev::gpd.fit"))
} else
if(package=="Renext") ##################
{
  if(is.null(method)) method <- 'r'
  if(method=="f")
  {
  outlist$q_gpd_creator <- "Renext::fGPD"
  z <- try(Renext::fGPD(x[x>0], ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "Renext::fGPD"))
  } else
  if(method=="r")
  {
  outlist$q_gpd_creator <- "Renext::Renouv"
  z <- try(Renext::Renouv(x, threshold=threshold, effDuration=length(x),
                    distname.y="gpd", plot=FALSE, ...), silent=TRUE)
  if(inherits(z, "try-error")) return(failfun(z, "Renext::Renouv"))
  } else
  stop("With package='Renext', method ('",method,"') must be 'f' or 'r'.")
} else
stop("package ", package, "is not in the options. This is a bug. Please report to berry-b@gmx.de.")
#
#
# quantile computing: ----------------------------------------------------------
if(package=="evir") ##################
{
  # computing part from evir::quant, Version: 1.7-3, Date: 2011-07-22
  lambda <- length(x)/z$n.exceed
  a <- lambda * (1 - probs)
  gfunc <- function(a, xihat) (a^(-xihat) - 1)/xihat
  output <- z$threshold + z$par.ests["beta"] * gfunc(a, z$par.ests["xi"])
} else
if(package=="evd") ##################
{
  probs2[probs2==0] <- NA
  probs2[probs2==1] <- NA
  output <- evd::qgpd(p=probs2, loc=z$threshold , scale=z$param["scale"], shape=z$param["shape"])
} else
if(package=="extRemes") ##################
{
  # Get parameters from result:
  if(z$method=="Bayesian")
  {
  scale <- exp(mean(z$results[,"log.scale"]))
  shape <-     mean(z$results[,"shape"])
  } else
  if(z$method=="Lmoments")
  {
  scale <- z$results["scale"]
  shape <- z$results["shape"]
  } else
  {
  scale <- z$results$par["scale"]
  shape <- z$results$par["shape"]
  }
  probs2[probs2==0] <- NA
  probs2[probs2==1] <- NA
  output <- extRemes::qevd(p=probs2, scale=scale, shape=shape, threshold=z$threshold, type="GP")
} else
if(package=="fExtremes") ##################
{
  output <- fExtremes::qgpd(p=probs2, xi=z@fit$par.ests["xi"], mu=z@parameter$u, beta=z@fit$par.ests["beta"])
  output <- as.vector(output)
  z2 <- lapply(slotNames(z), getElement, object=z)
  names(z2) <- slotNames(z)
  z <- z2
  outlist$q_gpd_Warning <- "transformed into list from Formal class 'fGPDFIT' [package 'fExtremes'] with 8 slots"
} else
if(package=="ismev") ##################
{
  # from ismev Version 1.40 Date 2009-14-07, Published 2014-12-24    ismev:::gpdq
  ismev_gpdq <- function(a,u,p) u + (a[1] * (p^(-a[2]) - 1))/a[2]
  output <- ismev_gpdq(a=z$mle, u=z$threshold, p=1-probs2)
} else
if(package=="Renext") ##################
{
  if(method=="f")
  output <- Renext::qGPD(probs2, scale=z$estimate["scale"], shape=z$estimate["shape"])
  else if(method=="r")
  output <- Renext::qGPD(probs2, loc=z$threshold, scale=z$estimate["scale"], shape=z$estimate["shape"])
} else
stop("package ", package, "is not in the options. This is a bug. Please report to berry-b@gmx.de.")
#
#
# output formatting, checks: ---------------------------------------------------
names(output) <- paste0(probs*100,"%")
# replace probs below truncation value with NA:
if(undertruncNA & any(probs < truncate) & !quiet)
  on.exit(message("Note in q_gpd: quantiles for probs (",
     toString(probs[probs<=truncate]),
     ") below truncate (",truncate,") replaced with NAs."), add=TRUE)
if(undertruncNA) output[probs < truncate] <- NA
# Output result:
outlist$q_gpd_method <- method
outlist$q_gpd_quant <- output
if(returnlist) output <- c(z, outlist)
output
}

