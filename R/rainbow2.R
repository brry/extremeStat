# reversed rainbow from blue to red
# Berry Boessenkool, Sept 2014

# calls \code{\link{rainbow} with different defaults and reverses the colors
rainbow2 <- function(
n=10, # number of colors
s=1, v=1, # saturation and value as in \code{\link{rainbow}}
start=0, # start color
end=0.7, # end color
alpha=1) # transparency
{
rev(rainbow(n=n, s=s, v=v, start=start, end=end))
}
