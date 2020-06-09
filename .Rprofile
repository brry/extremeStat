if(.Platform$OS.type=="unix") .libPaths("/home/berry/R/libBerry/")
if(requireNamespace("installB", quietly=TRUE))
{
if(.Platform$OS.type=="unix") installB::loadPackages(ask=FALSE)
installB::checkOutdated("extremeStat")
}
requireNamespace("testthat", quietly=TRUE)
utils::data(annMax)
