extremeStat
===========

Fit, plot and compare several (extreme value) distributions. 
Can also compute (truncated) distribution quantile estimates and draw a plot with return periods on a linear scale.

See the [Vignette](http://htmlpreview.github.io/?https://github.com/brry/extremeStat/blob/master/inst/doc/quantileEstimation.html) for an introduction to the package.

Install / load the package and browse through the examples:
```R
install.packages(c("lmomco", "pbapply","devtools"))
install.packages(c("evir", "ismev", "fExtremes", "extRemes", "evd", "Renext"))
# reiterate untill all of them work (some may not install properly on first try)

devtools::install_github("brry/berryFunctions")
devtools::install_github("brry/extremeStat") 

library(extremeStat)
?extremeStat
```
There are irregularly spaced github releases that follow the version number, eg
```R
devtools::install_github("brry/extremeStat@v0.4.38") # 2015-09
devtools::install_github("brry/extremeStat@v0.4.48") # 2016-02
```