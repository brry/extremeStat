# Function to get log-axis values and labels
# berry-b@gmx.de, Feb 2014, idea 2013

logVals <- function(
  from=-7,               # lower exponent OR vector with data
  to=7,                  # high end
  Range=range(from, to), # or give from and to as range
  base=1,                # bases to be used, eg. c(1,2,5)
  big.mark="'",          # symbol separating thousands, eg. space, comma, dot, etc. see "format" and "prettyNum"
  decimal.mark=".",      # character separating comma values, see "format" and "prettyNum"
  scientific=FALSE)      # see "format"
{
# Calculate the exponents from vector, if given as first argument:
if( missing(to)  &  NROW(from)>1  )
  {
  rng <- range(log10(from[from>0]), finite=TRUE)
  from <- floor(rng[1])
  to <- ceiling(rng[2])
  }
# or calculate the exponents from range, if given
if( !missing(Range)  )
  {
  from <- floor(Range[1])
  to <- ceiling(Range[2])
  }
# values for lines and labels:
vals <- base*10^rep(from:to, each=length(base))
# formatted values for labels:
labs <-  format(vals, big.mark=big.mark, trim=TRUE, scientific=scientific, 
                      drop0trailing=TRUE, decimal.mark=decimal.mark)
# Values for lines:
all <- 1:9 * 10^rep(from:to, each=9)
# return end result
list(vals=vals, labs=labs, all=all)
}
