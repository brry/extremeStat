### intro

Fit, plot and compare several (extreme value) distributions. 
Can also compute (truncated) distribution quantile estimates and draw a plot with return periods on a linear scale.

**See the [Vignette](https://cran.r-project.org/package=extremeStat/vignettes/extremeStat.html) for an introduction to the package.**

### installation

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/extremeStat)](https://cran.r-project.org/package=extremeStat) 
[![downloads](https://cranlogs.r-pkg.org/badges/extremeStat)](https://www.r-pkg.org/services)
[![Rdoc](https://www.rdocumentation.org/badges/version/extremeStat)](https://www.rdocumentation.org/packages/extremeStat)
!["extremeStat dependencies"](https://tinyverse.netlify.com/badge/extremeStat)


Install / load the package and browse through the examples:
```R
install.packages("extremeStat")
library(extremeStat)
vignette(extremeStat)

# update to the current development version, incl. vignette:
remotes::install_github("brry/extremeStat", build_vignettes=TRUE)
```

### trouble

If normal installation doesn't work, click on **Code - Download ZIP** (topright), unzip the file to some place, then
```R
setwd("that/path")
dd <- dir("extremeStat-master/R", full=T)
dummy <- sapply(dd, source)
```
This creates all R functions as objects in your globalenv workspace (and overwrites existing objects of the same name!).
