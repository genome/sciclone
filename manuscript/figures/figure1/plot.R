library(sciClone)
load("out.Rdata")
sc.plot1d(sc,"figure1.pdf", highlightSexChrs=TRUE, highlightsHaveNames=FALSE, overlayClusters=TRUE, showTitle=TRUE, cnToPlot=c(1:3))
