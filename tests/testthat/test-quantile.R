context("distLquantile")

data(annMax, package="extremeStat") # Annual Discharge Maxima (streamflow)
set.seed(007) # with other random samples, there can be warnings in q_gpd -> Renext::fGPD -> fmaxlo

ndist <- length(lmomco::dist.list()) - 13 + 22
# 13: excluded in distLfit.R Line 149
# 22: empirical, weighted, GPD_, n, threshold, etc

test_that("distLquantile generally runs fine",{
distLquantile(annMax)
expect_equal(nrow(distLquantile(annMax[annMax<30])), ndist)
expect_equal(nrow(distLquantile(annMax)), ndist)
expect_silent(distLquantile(annMax, truncate=0.6, gpd=FALSE, time=FALSE))
expect_message(distLquantile(annMax, selection="wak", empirical=FALSE, quiet=FALSE), 
  "distLfit execution took")
expect_message(distLquantile(rexp(199), truncate=0.8, probs=0.7, time=FALSE, emp=FALSE, quiet=FALSE),
  "must contain values that are larger than")
expect_message(distLquantile(rexp(4), selection="gpa"),
  "Note in distLquantile: sample size is too small to fit parameters (4). Returning NAs", fixed=TRUE)

d <- distLquantile(annMax, probs=0:4/4)
})


test_that("infinite values are removed",{
expect_message(distLextreme(c(-Inf,annMax)),
               "1 Inf/NA was omitted from 36 data points (2.8%)", fixed=TRUE)
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
expect_equal(nrow(distLquantile(annMax, emp=FALSE)), ndist-19) # only distributions in lmomco
aq <- distLquantile(annMax, truncate=0.8, probs=0.95) # POT
#round(aq,4)
# expected output (depending on lmomco version)
ex <- read.table(header=TRUE, text="
                           95%   RMSE
exp                   101.1631 0.0703
lap                   100.5542 0.0774
gpa                   103.4762 0.0778
wak                   103.4762 0.0778
wei                   102.7534 0.0796
pe3                   102.4791 0.0806
kap                   106.0260 0.0816
gno                   102.1442 0.0822
ln3                   102.1442 0.0822
gev                   101.9731 0.0831
glo                   101.4164 0.0870
pdq3                  101.2073 0.0875 # added Aug 2022
gum                   102.5499 0.0893
ray                   103.6840 0.0971
pdq4                  107.0252 0.1023 # added Aug 2022
gam                   103.8951 0.1128
rice                  104.2135 0.1217
nor                   104.2161 0.1218
revgum                104.9992 0.1595
empirical             109.2000     NA
quantileMean          105.7259     NA
weighted1             102.9910     NA # |
weighted2             102.8478     NA # | > changed Aug 2022, ignored in test
weighted3             102.5979     NA # | 
weightedc                  NaN     NA
GPD_LMO_lmomco        103.4762 0.0156
GPD_LMO_extRemes       99.8417 0.0163
GPD_PWM_evir          100.9874 0.0169
GPD_PWM_fExtremes     100.7009 0.0176
GPD_MLE_extRemes       99.0965 0.0161
GPD_MLE_ismev         108.8776 0.0467
GPD_MLE_evd           108.4444 0.0454
GPD_MLE_Renext_Renouv 108.4226 0.0453
GPD_MLE_evir                NA     NA
GPD_MLE_fExtremes           NA     NA
GPD_GML_extRemes      100.9103 0.0161 # changed from 99.0965 (2022-11-16) after bug fix by Eric G.
GPD_MLE_Renext_2par   166.9137 0.0958
GPD_BAY_extRemes            NA     NA
n_full                 35.0000     NA
n                       7.0000     NA
threshold              82.1469     NA")
colnames(ex) <- colnames(aq)
ex <- as.matrix(ex)

tsta <- rownames(aq) %in% lmomco::dist.list() | substr(rownames(aq),1,3) %in% c("GPD","n_f","n","thr")
tste <- rownames(ex) %in% lmomco::dist.list() | substr(rownames(ex),1,3) %in% c("GPD","n_f","n","thr")
tsta[rownames(aq)=="GPD_GML_extRemes"] <- FALSE # excluded while extRemes is being updated
tste[rownames(ex)=="GPD_GML_extRemes"] <- FALSE
if(is.na(aq["GPD_MLE_Renext_Renouv",1]))
{
tsta[rownames(aq)=="GPD_MLE_Renext_Renouv"] <- FALSE # excluded on weird Mac CRAN check
tste[rownames(ex)=="GPD_MLE_Renext_Renouv"] <- FALSE
}
expect_equal(round(aq[tsta,],1), round(ex[tste,],1))

dd <- distLquantile(annMax, selection="gpa", weighted=FALSE, truncate=0.001)
expect_equal(sum(is.na(dd[1:15,1:3])), 0)
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
               glob2rx("in pargno(lmom, ...): L-skew is too large*"), ignore.case=TRUE)

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


# wakeby (and others) with unrealistically high values:
xx6 <- c(0.342, 0.398, 0.415, 0.415, 0.462, 0.477, 0.491, 0.756, 0.763, 1.699)
d6 <- distLquantile(xx6, probs=c(0.8,0.9,0.99,0.9999), list=TRUE)
plotLfit(d6, xlim=c(0,2), nbest=10); d6$quant[1:10,]              # 36!!!
# works fine here:
xx7 <- c(0.415, 0.415, 0.431, 0.447, 0.531, 0.544, 0.643, 0.732, 0.82, 1.134)
d7 <- distLquantile(xx7, probs=c(0.8,0.9,0.99,0.9999), list=TRUE)
plotLfit(d7, xlim=c(0,2), nbest=10); d7$quant[1:10,]              # 4 (good)

})

