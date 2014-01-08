#!/usr/bin/env Rscript
library(sciClone)

tum = read.table("data/tumor.vafs", sep="\t", header=TRUE)
rel = read.table("data/relapse.vafs", sep="\t", header=TRUE)

highlight.genes <- read.table("data/genesToAnnotate", header=FALSE, sep="\t", as.is=TRUE)

samples = c("AML28 Tumor", "AML28 Relapse")

cn1 = read.table("data/tumor.cn", sep="\t")
cn1 = cn1[,c(1,2,3,5)]
cn2 = read.table("data/relapse.cn", sep="\t")
cn2 = cn2[,c(1,2,3,5)]

# Change this to 1 to output an intermediate plot at every iteration (as shown in supplement)
plotIntermediateResults <- 0

sc <- sciClone(vafs=list(tum,rel), sampleNames=samples, useSexChrs=FALSE, copyNumberCalls=list(cn1,cn2), copyNumberMargins=0.25, minimumDepth=100, verbose=1, doClusteringAlongMargins=TRUE, plotIntermediateResults = plotIntermediateResults)


writeClusterTable(sc, "clusters")
writeClusterSummaryTable(sc, "clusters.summary")
#sc.plot2d(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=positionsToHighlight, "figure3.pdf")
sc.plot2dWithMargins(sc, highlightsHaveNames=1, positionsToHighlight=highlight.genes, "figure3.pdf")

connectivity.matrix <- getConnectivityMatrix(sc)

write.table(file="connectivity.matrix.tsv", connectivity.matrix, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


save.image("out.Rdata")
