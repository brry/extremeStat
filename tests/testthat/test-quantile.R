context("distLquantile")

data(annMax, package="extremeStat") # Annual Discharge Maxima (streamflow)

test_that("distLquantile generally runs fine",{
distLquantile(annMax)
expect_equal(nrow(distLquantile(annMax)), 38)
expect_silent(distLquantile(annMax, truncate=0.6, gpd=FALSE, time=FALSE))
expect_message(distLquantile(annMax, selection="wak", empirical=FALSE, quiet=FALSE), 
  "distLfit execution took")
expect_message(distLquantile(rexp(199), truncate=0.8, probs=0.7, time=FALSE, emp=FALSE, quiet=FALSE),
  "must contain values that are larger than")
expect_message(distLquantile(rexp(4), selection="gpa"),
  "Note in distLquantile: sample size is too small to fit parameters (4). Returning NAs", fixed=TRUE)
})


test_that("distLquantile can handle selection input",{
dlf <- distLquantile(annMax, selection="wak", empirical=FALSE, list=TRUE)
plotLquantile(dlf, breaks=10)
expect_message(distLquantile(rexp(199), sel=c("wak", "gpa"), truncate=0.8, probs=c(0.7, 0.8, 0.9), quiet=FALSE),
  "Note in q_gpd: quantiles for probs (0.7, 0.8) below truncate (0.8) replaced with NAs.", fixed=TRUE)
distLquantile(rexp(199), selection=c("wak", "gpa"))
distLquantile(rexp(199), selection="gpa")
expect_error(distLquantile(rexp(199), selection=1:5, emp=FALSE), # index is a bad idea anyways
  "Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector", fixed=TRUE)
expect_error(distLquantile(rexp(199), selection=-3),
 "Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector", fixed=TRUE)

set.seed(42)
expect_warning(dlf <- distLfit(rnorm(100))) # gam + ln3 excluded
expect_equal(dlf$distfailed, c(gam="gam", ln3="ln3"))

dlf <- distLfit(annMax)
shouldbe <- c("80%"=82.002, "90%"=93.374, "99%"=122.505, "RMSE"=0.022)
expect_equal(distLquantile(annMax,  selection = "dummy"),
             distLquantile(dlf=dlf, selection = "dummy"))
d1 <- distLquantile(annMax,  selection = c("dummy","revgum","wak"))
d2 <- distLquantile(dlf=dlf, selection = c("dummy","revgum","wak"))
expect_equal(d1,d2)
expect_equal(round(d1[1,], 3), shouldbe)
expect_equal(round(d2[1,], 3), shouldbe)


dlf <- distLfit(annMax, selection=c("ln3","wak","gam", "gum"))
expect_equal(rownames(dlf$gof), c("wak", "ln3", "gum", "gam") )
sel <- c("dummy","gam","zzz","revgum","wak")
d3 <- distLquantile(annMax,  selection=sel, emp=FALSE )
d4 <- distLquantile(dlf=dlf, selection=sel, emp=FALSE )
o3 <- distLquantile(annMax,  selection=sel, emp=FALSE, order=FALSE)
o4 <- distLquantile(dlf=dlf, selection=sel, emp=FALSE, order=FALSE)

expect_equal(rownames(d3)[1:5], c("wak","gam","revgum","dummy","zzz"))
expect_equal(rownames(d4)[1:5], c("wak","gam","dummy","zzz","revgum")) # dlf does not have revgum
expect_equal(rownames(o3)[1:5], sel)
expect_equal(rownames(o4)[1:5], sel)
})


test_that("distLfit can handle truncate and threshold",{
  expect_message(dlf <- distLfit(annMax), "distLfit execution", all=TRUE)
  expect_message(dlf <- distLfit(annMax, truncate=0.7), "distLfit execution", all=TRUE)
  expect_message(dlf <- distLfit(annMax, threshold=50), "distLfit execution", all=TRUE)
  expect_message(dlf <- distLfit(annMax), "distLfit execution", all=TRUE)
})

test_that("distLquantile can deal with a given dlf",{
  dlf <- distLfit(annMax)
  expect_error(distLquantile(dlf, truncate=0.7), "x must be a vector")
  distLquantile(dlf=dlf, truncate=0.7)
  expect_message(dlf <- distLfit(annMax, threshold=50), "distLfit execution")
  expect_message(dlf <- distLfit(annMax), "distLfit execution")
})

test_that("distLquantile can handle emp, truncate",{
expect_equal(nrow(distLquantile(annMax, emp=FALSE)), 20) # only distributions in lmomco
aq <- distLquantile(annMax, truncate=0.8, probs=0.95) # POT
expect_equal(round(mean(aq[1:35,1], na.rm=TRUE),2), 102.1) 
dd <- distLquantile(annMax, selection="gpa", weighted=FALSE, truncate=0.001)
expect_equal(sum(is.na(dd[,1:3])), 3)
expect_equal(dd["gpa",1:3], dd["GPD_LMO_lmomco",1:3])
})

test_that("distLquantile can handle list",{
# Compare several GPD Fitting functions:
distLquantile(annMax, threshold=70, selection="gpa", weighted=FALSE, list=TRUE)
expect_is(distLquantile(annMax, truncate=0.62, list=TRUE), "list")
expect_is(distLquantile(annMax, threshold=70,  list=TRUE), "list")
})


