context("distLquantile")

data(annMax, package="extremeStat") # Annual Discharge Maxima (streamflow)


test_that("distLquantile generally runs fine",{
distLquantile(annMax)
expect_silent(distLquantile(annMax, truncate=0.6, gpd=FALSE, quiet=TRUE))
expect_message(distLquantile(annMax, plot=FALSE, selection="wak", empirical=FALSE), 
  "Note in distLgof: Only wak was fitted, thus GOF can't be compared.")
distLquantile(annMax, plot=TRUE, selection="wak", empirical=FALSE, breaks=10)
expect_message(distLquantile(rexp(199), sel=c("wak", "gpa"), truncate=0.8, probs=c(0.7, 0.8, 0.9)),
  "Note in q_gpd: quantiles for probs (0.7, 0.8) below truncate (0.8) replaced with NAs.", fixed=TRUE)
expect_message(distLquantile(rexp(199), truncate=0.8, probs=0.7, time=FALSE, emp=FALSE),
  "must contain values that are larger than")
distLquantile(rexp(199), selection=c("wak", "gpa"))
distLquantile(rexp(199), selection="gpa")
distLquantile(rexp(4))
expect_warning(distLquantile(rexp(4), selection="gpa"),
               "gpa available in dlf$parameter, but not in dlf$gof", fixed=TRUE)
expect_error(distLquantile(rexp(199), selection=1:5, emp=FALSE), # index is a bad idea anyways
  "Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector", fixed=TRUE)
expect_error(distLquantile(rexp(199), selection=-3),
 "Since Version 0.4.36 (2015-08-31), 'selection' _must_ be a character string vector", fixed=TRUE)
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
expect_equal(nrow(distLquantile(annMax, emp=FALSE)), 17) # only distributions in lmomco
expect_message(aq <- distLquantile(annMax, truncate=0.8, probs=0.95)) # POT
expect_equal(mean(aq, na.rm=TRUE), 102.055973)
expect_equal(sum(is.na(distLquantile(annMax, selection="gpa", weight=FALSE))), 3)
})

test_that("distLquantile can handle returndlf",{
# Compare several GPD Fitting functions:
distLquantile(annMax, threshold=70, selection="gpa", weight=FALSE, returndlf=TRUE)
expect_is(distLquantile(annMax, truncate=0.62, returndlf=TRUE), "list")
expect_is(distLquantile(annMax, threshold=70,  returndlf=TRUE), "list")
})


