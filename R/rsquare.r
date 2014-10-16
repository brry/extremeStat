rsquare <- function(
                    a,
                    b,
                    quiet=FALSE)
{
if(!(is.vector(a) & is.vector(b))) stop("input is not vectors")
if(length(a) != length(b)) stop("vectors not of equal length")
if(any(is.na(a)|is.na(b)))
     {
     Na <- which(is.na(a)|is.na(b))
     if(!quiet) warning(length(Na), " NAs were omitted from ", length(a), " data points.")
     a <- a[-Na] ; b <- b[-Na]
     } # end if NA
cor(a,b)^2
}



if(FALSE)
{
# alternative, slower (3.4 instead of 2.1 seconds in the example below)
# crucial, if calculations are done iteratively or performed multiple times
rsquare2 <- function(a,b) { 
  if(!(is.vector(a) & is.vector(b))) stop("input is not vectors")
  if(length(a) != length(b)) stop("vectors not of equal length")
  if(any(is.na(a)|is.na(b)))
     { warning("NAs were omitted")
     Na <- which(is.na(a)|is.na(b))
     a <- a[-Na] ; b <- b[-Na]
     } # end if NA
  aa <-  a-mean(a)
  bb <-  b-mean(b)
  sum(aa*bb)^2/sum(aa^2)/sum(bb^2) }

a <- sort(rnorm(1e8)); b <- 2*a+3+rnorm(length(a))
system.time(rsquare(a,b))
system.time(rsquare2(a,b))
}
