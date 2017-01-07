context("distLquantile")

data(annMax, package="extremeStat") # Annual Discharge Maxima (streamflow)
set.seed(007) # with other random samples, there can be warnings in q_gpd -> Renext::fGPD -> fmaxlo


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

d <- distLquantile(annMax, probs=0:4/4)
})


test_that("distLquantile can handle selection input",{
dlf <- distLquantile(annMax, selection="wak", empirical=FALSE, list=TRUE)
plotLquantile(dlf, breaks=10)
expect_message(distLquantile(rexp(199), sel=c("wak", "gpa"), truncate=0.8, probs=c(0.7, 0.8, 0.9)),
  "Note in q_gpd: quantiles for probs (0.7) below truncate (0.8) replaced with NAs.", fixed=TRUE)
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
d1 <- distLquantile(annMax,  selection="dummy", onlydn=FALSE)
d2 <- distLquantile(dlf=dlf, selection="dummy", onlydn=FALSE)
expect_equal(d1,d2)
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
expect_equal(round(aq[1:35,1],2), 
             structure(c(101.16, 100.55, 103.48, 103.48, 102.75, 102.48, 106.03, 
102.14, 102.14, 101.97, 101.42, 102.55, 103.68, 103.9, 104.21, 
104.22, 105, 105.73, 102.88, 102.74, 102.66, NaN, 103.48, 99.84, 
100.99, 100.7, 99.1, 108.88, 108.44, 108.42, NA, NA, 99.1, 166.91,  # Renext2par 80.25
NA), .Names = c("exp", "lap", "gpa", "wak", "wei", "pe3", "kap", 
"gno", "ln3", "gev", "glo", "gum", "ray", "gam", "rice", "nor", 
"revgum", "quantileMean", "weighted1", "weighted2", "weighted3", 
"weightedc", "GPD_LMO_lmomco", "GPD_LMO_extRemes", "GPD_PWM_evir", 
"GPD_PWM_fExtremes", "GPD_MLE_extRemes", "GPD_MLE_ismev", "GPD_MLE_evd", 
"GPD_MLE_Renext_Renouv", "GPD_MLE_evir", "GPD_MLE_fExtremes", 
"GPD_GML_extRemes", "GPD_MLE_Renext_2par", "GPD_BAY_extRemes"
))) 
dd <- distLquantile(annMax, selection="gpa", weighted=FALSE, truncate=0.001)
expect_equal(sum(is.na(dd[1:15,1:3])), 3)
expect_equal(dd["gpa",1:3], dd["GPD_LMO_lmomco",1:3])
})

test_that("distLquantile can handle list",{
# Compare several GPD Fitting functions:
distLquantile(annMax, threshold=70, selection="gpa", weighted=FALSE, list=TRUE)
expect_is(distLquantile(annMax, truncate=0.62, list=TRUE), "list")
expect_is(distLquantile(annMax, threshold=70,  list=TRUE), "list")
})


test_that("distLquantile can handle inputs with (rare) errors",{
# invalid lmoms
xx1 <- c(4.2, 1.1, 0.9, 5, 0.6, 5.1, 0.9, 1.2, 0.6, 0.7, 0.9, 1.1, 1.3, 
1.4, 1.4, 0.6, 3, 1.6, 0.5, 1.4, 1.1, 0.5, 1.3, 3.6, 0.5)
expect_message(distLquantile(xx1, truncate=0.8), 
               "Note in distLfit: L-moments are not valid. No distributions are fitted.")

# kap failed
xx2 <- c(0.6, 1.6, 2.2, 0.6, 0.9, 3.3, 1.3, 4.7, 0.9, 0.8, 0.5, 0.8, 0.6, 0.7, 1.1, 0.9, 
       5.4, 3.9, 0.9, 0.7, 0.6, 0.7, 15.1, 2.7, 0.7, 1, 0.5, 0.6, 1, 0.9, 1.4)
dd <- distLquantile(xx2, truncate=0.8)
expect_equal(dd["kap","RMSE"], NA_real_)

# kap and ln3
xx3 <- c(0.7, 1.5, 0.7, 2.6, 0.7, 0.8, 1.9, 5.4, 1.4, 1, 1.7, 0.8, 1.3, 0.8, 0.9, 0.5, 
       0.5, 5.1, 0.9, 1, 1, 1.4, 1.5, 1.4, 4.9, 0.6, 4.3, 0.7, 0.7, 1.2, 0.9, 0.8)
expect_warning(dd <- distLquantile(xx3, truncate=0.8), 
               glob2rx("in parln3(lmom, ...): L-skew is negative, try reversing the data*"))
expect_equal(dd["kap","RMSE"], NA_real_)

# strongly skewed (gno):
xx4 <- c(2.4,2.7,2.3,2.5,2.2, 62.4 ,3.8,3.1) 
expect_warning(dd <- distLquantile(xx4), 
               glob2rx("in pargno(lmom, ...): L-SKEW IS TOO LARGE FOR ROUTINE*"))

# kap should fail:
xx5 <- c(2.4, 2.5, 2.6, 2.9, 4.2, 4.6, 5.7)
distLfit(xx5)$parameter$kap
dfun <- function(xxx) expect_true(all(is.na(distLquantile(xxx, probs=0:10/10, 
                                            sel="kap", emp=FALSE)["kap",])))
dfun(xx5)
dfun(c(2.2, 2.3, 2.3, 2.3, 4.1, 8.8))
dfun(c(2.2, 2.3, 2.4, 2.5, 3.2, 4.2, 4.5, 5.9, 6))
dfun(c(1.8, 1.8, 2, 2, 2.6, 2.7, 3.7, 3.7))
dfun(c(2.2, 2.2, 2.3, 2.9, 3.4, 4.4, 5.2))
dfun(c(2.1, 2.2, 2.5, 3.2, 7.8, 16.1)) # kap has 4 distinct values here...
})
