library(sciClone)
load("out.Rdata")


sc.plot1d(sc,"figure2a.pdf", highlightSexChrs=FALSE, highlightsHaveNames=FALSE, overlayClusters=TRUE, showTitle=FALSE, cnToPlot=c(2), biggerText=TRUE)

sc.plot1d(sc2,"figure2b.pdf", highlightSexChrs=FALSE, highlightsHaveNames=FALSE, overlayClusters=TRUE, showTitle=FALSE, cnToPlot=c(2), biggerText=TRUE)

