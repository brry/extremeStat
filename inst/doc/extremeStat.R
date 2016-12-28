## ----instcran, eval=FALSE------------------------------------------------
#  install.packages("extremeStat")
#  library(extremeStat)

## ----instgit, eval=FALSE-------------------------------------------------
#  install.packages(c("devtools","evd","evir","extRemes","fExtremes",
#                     "ismev","lmomco","pbapply","Renext"))
#  # reiterate untill all of them work (some may not install properly on first try)
#  
#  devtools::install_github("brry/berryFunctions")
#  devtools::install_github("brry/extremeStat")
#  library(extremeStat)

## ----library, echo=FALSE-------------------------------------------------
library(extremeStat)

## ----dataHist, fig.show='hold', echo=-3----------------------------------
data(rain, package="ismev")
rain <- rain[rain>2]
par(mar=c(3.2,3.2,1.5,0.7), mgp=c(2.1,0.7,0))
hist(rain, breaks=80, col=4, las=1)
# Visual inspection is easier on a logarithmic scale:
berryFunctions::logHist(rain, breaks=80, col=3, las=1)

## ----dlq-----------------------------------------------------------------
dlq <- distLquantile(rain[1:900], probs=c(0.8,0.9,0.99,0.999), list=TRUE, quiet=TRUE)

## ----dlprint, eval=1-----------------------------------------------------
printL(dlq)
# More information on dlf objects in
?extremeStat

## ----dlplot, echo=-1, fig.height=3.5, fig.width=5.5----------------------
par(mar=c(3.9,3.9,1.5,0.7), mgp=c(2.8,0.7,0))
plotLquantile(dlq, nbest=8, linargs=list(lwd=2), 
          heights=seq(0.04, 0.01, len=8), breaks=80)

## ----dlquant-------------------------------------------------------------
dlq$quant # distLquantile output if returnlist=FALSE (the default)

## ----weight, echo=-1, fig.height=3.5, fig.width=5.5----------------------
par(mar=c(3.2,3.6,2.6,0.7), mgp=c(2.1,0.7,0))
plotLgof(dlq, legargs=list(cex=0.8, bg="transparent") )

## ----prob, echo=FALSE, fig.height=1, fig.width=5.5-----------------------
  par(mar=rep(0,4))
  plot(1:6, type="n", ylim=c(1,6), axes=F, ann=F)
  arrows(x0=1, y0=4, x1=6, code=3, angle=90, length=0.07)
  arrows(x0=4, y0=4, x1=6, code=1, angle=90, length=0.07)
  arrows(x0=4, y0=3, x1=6, code=3, angle=90, length=0.07)
  text(x=c(1,4,5.5,6), y=rep(5,4), c(0,"t","p",1), col=c(1,1,2,1))
  text(x=c(4,5.5,6),   y=rep(2,3), c(0,"p2",1),    col=c(1,2,1))
  text(3.8, 3, "truncated sample", adj=1)
  points(x=rep(5.5,2), y=c(3,4), col=2, pch=16)
#text(7, 5.5, expression(frac(1-p, 1-t)*"  =  "* frac(1-p2, 1-0)), adj=0, cex=1.2)
#text(7, 2.5, expression("p2  =  "*frac(p-t, 1-t)), adj=0, cex=1.2)


## ----trunc, echo=-1, fig.height=3.5, fig.width=5.5-----------------------
par(mar=c(3.2,3.6,2.6,0.7), mgp=c(2.1,0.7,0))
d <- distLquantile(rain, truncate=0.9, plot=TRUE, probs=0.999, quiet=TRUE, breaks=50)

## ----teff, eval=FALSE----------------------------------------------------
#  tt <- seq(0,0.95, len=50)
#  if(interactive()) lapply <- pbapply::pblapply # for progress bars
#  qq <- lapply(tt, function(t) distLquantile(rain, truncate=t,
#                                               probs=c(0.99,0.999), quiet=TRUE))
#  save(tt,qq, file="qq.Rdata")

## ----teffplot, fig.height=3.5, fig.width=5.5-----------------------------
load("qq.Rdata")
par(mar=c(3,2.8,2.2,0.4), mgp=c(1.8,0.5,0))
plot(tt,tt, type="n", xlab="truncation proportion", ylab="Quantile estimation",
     main="truncation effect for 6k values of rain", ylim=c(22,90), las=1)
dn <- c("wak","kap","wei","gpa","pe3","weighted2")
cols <- c(4,5,3,"orange",2,1) ; names(cols) <- dn
for(d in rownames(qq[[1]])) lines(tt, sapply(qq, "[", d, j=2), col=8)
for(d in dn)
  {
  lines(tt, sapply(qq, "[", d, j=1), col=cols[d], lwd=2)
  lines(tt, sapply(qq, "[", d, j=2), col=cols[d], lwd=2)
  }
abline(h=berryFunctions::quantileMean(rain, probs=c(0.99,0.999)), lty=3)
legend("topright", c(dn,"other"), col=c(cols,8), lty=1, lwd=c(rep(2,6),1), bg="white", cex=0.6)
text(0.9, 53, "Q99.9%") ; text(0.9, 34, "Q99%")
text(0.35, 62, "empirical quantile (full sample)", cex=0.7)


## ----ssdep, eval=FALSE---------------------------------------------------
#  set.seed(1)
#  ss <- c(30,50,70,100,200,300,400,500,1000)
#  rainsamplequantile <- function() sapply(ss, function(s) distLquantile(sample(rain,s),
#            probs=0.999, plot=F, truncate=0.8, quiet=T, sel="wak", gpd=F, weight=F))
#  sq <- pbapply::pbreplicate(n=100, rainsamplequantile())
#  save(ss,sq, file="sq.Rdata")

## ----ssdepplot, fig.height=3.5, fig.width=5.5----------------------------
load("sq.Rdata")
par(mar=c(3,2.8,2.2,0.4), mgp=c(1.7,0.5,0))
sqs <- function(prob,row) apply(sq, 1:2, quantile, na.rm=TRUE, probs=prob)[row,]
berryFunctions::ciBand(yu=sqs(0.6,1), yl=sqs(0.4,1), ym=sqs(0.5,1), x=ss, 
    ylim=c(25,75), xlim=c(30,900), xlab="sample size", ylab="estimated 99.9% quantile", 
    main="quantile estimations of small random samples", colm="blue")
berryFunctions::ciBand(yu=sqs(0.6,2), yl=sqs(0.4,2), ym=sqs(0.5,2), x=ss, add=TRUE)
abline(h=quantile(rain,0.999))
text(250, 50, "empirical", col="forestgreen")
text(400, 62, "Wakeby", col="blue")
text(0, 61, "'True' population value", adj=0)
text(600, 40, "median and central 20% of 100 simulations")

## ----RP, echo=-1, warning=FALSE, fig.height=4, fig.width=5.5-------------
par(mar=c(3,2.8,1.2,0.4), mgp=c(1.8,0.5,0))
data("annMax") # annual discharge maxima in the extremeStat package itself
dle <- distLextreme(annMax, log=TRUE, legargs=list(cex=0.6, bg="transparent"), nbest=17, quiet=TRUE)
dle$returnlev[1:20,]


## ----help, eval=FALSE----------------------------------------------------
#  ?extremeStat

