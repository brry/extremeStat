lim0 <- function(
                 x,
                 f=1/27,
                 curtail=TRUE)
     # returns vector with 2 values: 0 and by 4% contra-extended max (as for xaxs="r")
     # see methods(plot), graphics:::plot.default,
     #     graphics:::curve and extendrange in grDevices
     {
     r <- range(x, finite=TRUE)
     r2 <- r + c(-f,f) * diff(r) # classical procedure of extendrange
     r2[which.min(abs(r2))] <- 0 # set one end to zero
     if(curtail) # if par xaxs is "r" as it is by default, first trim the range, so that
     r2 + c(f,-f) * diff(r2) # in the plot command, we don't have to change yaxs or xaxs
     else
     r2
     }
