#!/usr/bin/env Rscript

library(sciClone)


#BRC sample
rcnt = read.table(paste("data/TCGA-A2-A0YG/vafs",sep=""),sep="\t")
cn = read.table(paste("data/TCGA-A2-A0YG/cn",sep=""),sep="\t")

sc = sciClone(vafs=list(rcnt), sampleNames=c("TCGA−A2−A0YG"), copyNumberCalls=list(cn),
  cnCallsAreLog2=TRUE, minimumDepth=50)
writeClusterTable(sc,"clusters.brc")
writeClusterSummaryTable(sc,"clusters.summary.brc")

#UCEC sample
rcnt = read.table("data/TCGA-D1-A17T/variants.rcnt")
cn0 = read.table("data/TCGA-D1-A17T/copy_number.cbs")
cn0 = cn0[,c(1,2,3,5)]

sc2 = sciClone(vafs=list(rcnt), sampleNames=c("TCGA-D1-A17T"), copyNumberCalls=list(cn0), minimumDepth=50, cnCallsAreLog2=TRUE)
writeClusterTable(sc2, "clusters.ucec")
writeClusterSummaryTable(sc2, "cluster.summary.ucec")



#output cluster probabilities for genes of interest
anno = read.table("data/TCGA-D1-A17T/genesToAnnotate")
b = read.table("clusters.ucec",header=T)

cat("gene\tcluster\tcluster.prob.1\tcluster.prob.2\n")
for(i in 1:length(anno[,1])){
  x=b[b$chr==anno[i,1] & b$st==anno[i,2],];
  cat(paste(anno[i,3],x$cluster,x$cluster.prob.1,x$cluster.prob.2,sep="\t"),"\n")
}


save.image("out.Rdata")

source("plot.R")
