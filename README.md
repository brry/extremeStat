extremeStat
===========

Fit, plot and compare several (extreme value) distributions by means of a plot with return periods on a linear scale

independent from berryFunctions now. Functions that currently exist twice:
lim0, logAxis, logVals, owa, rainbow2, rmse, rsquare

In my opinion, it would make more sense to depend on berryFunctions, but Ripley does not like to see a subpackage like this on CRAN.
What's your opinion?


if(!require(devtools)) {install.packages("devtools"); require(devtools)}

install_github("BerryBoessenkool/extremeStat")

library(extremeStat)

?extremeStat

