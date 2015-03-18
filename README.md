extremeStat
===========

Fit, plot and compare several (extreme value) distributions and draw a plot with return periods on a linear scale.

depends on berryFunctions for:
lim0, logAxis, logVals, owa, rainbow2, rmse, rsquare

Ripley does not like to see a subpackage like this on CRAN.
What's your opinion?

Code to install:

```R
if(!require(devtools)) install.packages("devtools")
devtools::install_github("BerryBoessenkool/berryFunctions")
devtools::install_github("BerryBoessenkool/extremeStat")
library(extremeStat)
?extremeStat
```
