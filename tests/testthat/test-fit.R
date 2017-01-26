context("Fit distributions")

suppressPackageStartupMessages(library(lmomco))
data("annMax")

test_that("plotLfit ordering works fine",{
dlf <- distLfit(annMax)
sel <- c(pe3="pe3", revgum="revgum", rice="rice")

d <- plotLfit(dlf, cdf=TRUE, sel=sel)
printL(d)
expect_equal(d$distnames, sel[c(1,3,2)])

d <- plotLfit(dlf, cdf=TRUE, sel=sel, order=F)
expect_equal(d$distnames, sel)
})


if(FALSE){ 

data(rain, package="ismev")
samp <- function(n, seed=42) {set.seed(seed); sample(rain[rain>0.1], n)}
dlf <- distLfit(samp(25))
plotLfit(dlf, cdf=T)
plotLfit(  distLfit(samp(25,1))  )

res <- pblapply(1:1000, function(s) tryStack( distLfit(samp(25,s), quiet=TRUE) )  )

failed <- sapply(res, "[[", "distfailed")
table(failed)

}
