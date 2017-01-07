context("Goodness of fit + distribution weights")

test_that("distLweights ordering works fine",{
rn <- function(...) rownames(distLweights(...))
vect <- c(gum=0.20, wak=0.17, gam=0.21)
expect_equal(rn(vect),c("wak","gum","gam"))
expect_equal(rn(vect, order=FALSE),c("gum", "wak", "gam"))
vect[3] <- NA
expect_equal(rn(vect),c("wak","gum","gam"))
expect_equal(rn(vect, order=FALSE),c("gum", "wak", "gam"))
set.seed(007)
expect_equal(distLweights(data.frame(RMSE=sample(1:12)),onlydn=FALSE)[7:12,4], rep(0,6))
expect_equal(distLweights(data.frame(RMSE=sample(1:12)),onlydn=FALSE, order=FALSE)[1,4], 0)
})

test_that("distLweights handles single values",{
expect_equal(distLweights(c(gum=0.20)), data.frame(RMSE=0.2,
                                                   weight1=1,
                                                   weight2=1,
                                                   weight3=1,
                                                   weightc=0, row.names="gum"))
})

test_that("distLweights handles NAs",{
expect_equal(distLweights(c(gum=0.20, gam=NA, gev=NA))[,3], c(1,0,0))
expect_equal(distLweights(c(gum=0.20, gam=NA, gev=NA))[,4], c(1,0,0))
expect_equal(distLweights(c(gum=0.20, gam=NA, gev= 1))[,4], c(1,0,0))
})

test_that("distLweights can get RMSE from data.frame",{
  
set.seed(42); x <- data.frame(A=1:5, RMSE=runif(5)) ; x
expect_warning(distLweights(x), "There are no distributions matching lmomco::dist.list()")

expect_equal(round(distLweights(x,onlydn=FALSE),3),
structure(list(RMSE = c(0.286, 0.642, 0.83, 0.915, 0.937), weight1 = c(0.374, 
0.232, 0.157, 0.123, 0.114), weight2 = c(0.605, 0.275, 0.099, 
0.021, 0), weight3 = c(0.618, 0.28, 0.101, 0, 0), weightc = c(0, 
0, 0, 0, 0)), .Names = c("RMSE", "weight1", "weight2", "weight3", 
"weightc"), row.names = c("3", "5", "4", "1", "2"), class = "data.frame"))

expect_equal(nrow(distLweights(data.frame(RMSE=1:2),onlydn=FALSE)),2)
expect_equal(nrow(distLweights(data.frame(Rmse=1:3),onlydn=FALSE)),3)
expect_equal(nrow(distLweights(data.frame(rmse=1:4),onlydn=FALSE)),4)
expect_error(distLweights(data.frame(RMS=1:4)), "There is no column")
expect_error(distLweights(data.frame(A=1:4, RMSE=1:4, rmse=2:5)), "There are several columns")
})

test_that("distLweights RMSE must have names",{
expect_error(distLweights(1:8), "RMSE must have names")
})

test_that("distLweights weight3 has half of values zero",{
nzeros <- function(n, ...) sum(distLweights(data.frame(RMSE=1:n), onlydn=FALSE, ...)[,"weight3"]==0)
expect_equal(nzeros(1), 0)
expect_equal(nzeros(2), 1)
expect_equal(nzeros(3), 1)
expect_equal(nzeros(4), 2)
expect_equal(nzeros(6), 3)
expect_equal(nzeros(11), 5)
expect_equal(nzeros(12), 6)
})


test_that("distLweights checks for unused arguments",{
expect_warning(distLquantile(rexp(199), ks=FALSE, hi=TRUE), "unused arguments in")
})

