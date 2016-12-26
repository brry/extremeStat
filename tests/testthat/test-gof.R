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
expect_equal(distLweights(data.frame(rmse=sample(1:12)))[7:12,4], rep(0,6))
expect_equal(distLweights(data.frame(rmse=sample(1:12)), order=FALSE)[1,4], 0)
})

test_that("distLweights handles single values",{
expect_equal(distLweights(c(gum=0.20)), data.frame(rmse=0.2,
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
expect_equal(nrow(distLweights(data.frame(rmse=1:2))),2)
expect_equal(nrow(distLweights(data.frame(Rmse=1:3))),3)
expect_equal(nrow(distLweights(data.frame(RMSE=1:4))),4)
expect_error(distLweights(data.frame(RMS=1:4)), "There is no column")
expect_error(distLweights(data.frame(A=1:4, RMSE=1:4, rmse=2:5)), "There are several columns")
})

test_that("distLweights rmse must have names",{
expect_error(distLweights(1:8), "rmse must have names")
})

test_that("distLweights weight3 has half of values zero",{
nzeros <- function(n, ...) sum(distLweights(data.frame(rmse=1:n), ...)[,"weight3"]==0)
expect_equal(nzeros(1), 0)
expect_equal(nzeros(2), 1)
expect_equal(nzeros(3), 1)
expect_equal(nzeros(4), 2)
expect_equal(nzeros(6), 3)
expect_equal(nzeros(11), 5)
expect_equal(nzeros(12), 6)
})
