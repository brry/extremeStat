extremeStat
===========

Fit, plot and compare several (extreme value) distributions. 
Can also compute (truncated) distribution quantile estimates and draw a plot with return periods on a linear scale.


**See the [Vignette](https://cran.r-project.org/web/packages/extremeStat/vignettes/extremeStat.html) for an introduction to the package.**


[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/extremeStat)](http://cran.r-project.org/package=extremeStat) [![downloads](http://cranlogs.r-pkg.org/badges/extremeStat)](http://www.r-pkg.org/services)

Install / load the package and browse through the examples:
```R
install.packages("extremeStat")
library(extremeStat)
vignette(extremeStat)
```
There are irregularly spaced github releases that follow the version number:
```R
devtools::install_github("brry/extremeStat@v0.4.38") # 2015-09
devtools::install_github("brry/extremeStat@v0.4.48") # 2016-02
devtools::install_github("brry/extremeStat@v0.5.13") # 2016-03
```

Code to install the most recent development version from github:

```R
# Avoid installing devtools with all its dependencies:
source("http://raw.githubusercontent.com/brry/berryFunctions/master/R/instGit.R")
instGit("brry/extremeStat")

# or using devtools:
if(!requireNamespace("devtools", quitly=TRUE)) install.packages("devtools")
devtools::install_github("brry/extremeStat")

library(extremeStat)
```

If direct installation from CRAN doesn't work, your R version might be too old. In that case, an update is really recommendable: [r-project.org](http://www.r-project.org/). If you can't update R, try installing from source (github) via `instGit` or devtools as mentioned above. If that's not possible either, here's a manual workaround:
click on **Download ZIP** (to the right, [link](https://github.com/brry/extremeStat/archive/master.zip)), unzip the file to some place, then
```R
setwd("that/path")
dd <- dir("extremeStat-master/R", full=T)
dummy <- sapply(dd, source)
```
This creates all R functions as objects in your globalenv workspace (and overwrites existing objects of the same name!).
