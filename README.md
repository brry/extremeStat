extremeStat
===========

Fit, plot and compare several (extreme value) distributions. 
Can also compute (truncated) distribution quantile estimates and draw a plot with return periods on a linear scale.

Code to install:

```R
if(!require(devtools)) install.packages("devtools")
devtools::install_github("brry/berryFunctions")
devtools::install_github("brry/extremeStat")
library(extremeStat)
?extremeStat
```
There are irregularly spaced github releases that follow the version number, eg
```R
devtools::install_github("brry/extremeStat@v0.4.38")
```