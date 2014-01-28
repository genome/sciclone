#!/usr/bin/env Rscript

library(sciClone)

v0 = read.table("data/mmy.snv.vafs", sep="\t")
cn0 = read.table("data/cn.dat")
cn0 = cn0[,c(1,2,3,5)]
reg0 = read.table("data/exclude.loh")
reg0 = reg0[,c(1,2,3)]
clusterParams="empty"
#clusterParams="no.pvalue.outlier.detection"
sc = sciClone(vafs=list(v0), sampleNames=c("MMY4"), copyNumberCalls=list(cn0), regionsToExclude=list(reg0), minimumDepth=100, doClustering=TRUE, clusterParams=clusterParams,maximumClusters=10, copyNumberMargins=0.25, useSexChrs=FALSE)
writeClusterTable(sc, "clusters")
writeClusterSummaryTable(sc, "cluster.summary")
save.image("out.Rdata")
source("plot.R")
