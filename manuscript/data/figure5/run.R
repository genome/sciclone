#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(sciClone))

#read in data
regions_to_exclude = read.table(file='data/loh.regions',sep="\t",header=F);

samples = c("Pre-treatment Tumor 1","Pre-treatment Tumor 2","Post-treatment Tumor")

v1 = read.table("data/pre1.vafs",sep="\t")
v2 = read.table("data/pre2.vafs",sep="\t")
v3 = read.table("data/post1.vafs",sep="\t")

cn1 = read.table("data/pre1.cn",sep="\t")
cn1 = cn1[,c(1,2,3,5)]
cn2 = read.table("data/pre2.cn",sep="\t")
cn2 = cn2[,c(1,2,3,5)]
cn3 = read.table("data/post1.cn",sep="\t")
cn3 = cn3[,c(1,2,3,5)]

#cluster
sc = sciClone(vafs=list(v1,v2,v3), sampleNames=samples, copyNumberCalls=list(cn1,cn2,cn3), doClusteringAlongMargins=FALSE, maximumClusters=10, regionsToExclude=regions_to_exclude)
writeClusterTable(sc,"cluster")
writeClusterSummaryTable(sc,"cluster.summary")

save.image("out.Rdata")
source("plot.R")
