#!/usr/bin/env Rscript

# This file creates:
# Fig4A:  the grand consensus matrix (over binomial, beta, and gaussian, as
#         we vary number of clusters)
# Fig4B:  simulation results showing the number of clusters as we vary the
#         number of data points.
# Supp2:  the binomial 2D clustering of a characteristic clustering that differs
#         from the beta clustering shown in Fig 3 (initialized with 10 clusters).
# Supp1:  the gaussian 2D clustering of a characteristic clustering that differs
#         from the beta clustering shown in Fig 3 (initialized with 10 clusters).

suppressPackageStartupMessages(library(sciClone))

#library(stats)
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))

source("functions.R")

#col.palette <- brewer.pal(11,"RdBu")
#col.palette <- bluered(5)

v1 = read.table("data/tumor.vafs", sep="\t", header=TRUE)
v2 = read.table("data/relapse.vafs", sep="\t", header=TRUE)

genesToAnnotate <- read.table("data/genesToAnnotate", header=FALSE, sep="\t", as.is=TRUE)
highlightsHaveNames=TRUE

samples = c("AML28 Tumor", "AML28 Relapse")

cn1 = read.table("data/tumor.cn", sep="\t")
cn1 = cn1[,c(1,2,3,5)]
cn2 = read.table("data/relapse.cn", sep="\t")
cn2 = cn2[,c(1,2,3,5)]


# Use the results from fig 3 as the gold standard.

# Read in the clustering results
cluster.results <- read.table(file="data/clusters", sep="\t", header=TRUE)
num.clusters.orig <- length(unique(cluster.results$cluster))

# Determine the points that have adequate depth and are in copy number
# neutral regions.  And were not deemed outliers (i.e., in cluster 0).

valid.rows <- apply(cluster.results[grep("^cleancn",names(cluster.results))], 1, function(row) all(row==2)) & cluster.results$adequateDepth == 1 & !is.na(cluster.results$cluster & (cluster.results$cluster != 0))

valid.results <- cluster.results[valid.rows,]

valid.results.chr.start <- valid.results[,c(1,2)]






# Create Fig3A--the binomial consensus matrix as the number of clusters
# are varied.

# Some margin adjustments for the figures
#lhei <- c(0.08,0.92)
#lwid <- c(0.05,0.95)
#lhei <- c(0.01,0.92)
#lwid <- c(0.01,0.95)
#margins <- c(0.01,0.01)
#lhei <- c(0.00001,0.92)
#lwid <- c(0.00001,0.95)
#margins <- c(0,0)

# Put legend above heatmap
lmat <- rbind(c(0,4),c(2,1),c(0,3))
lhei <- c(1.0, 4.5, 0.1)
lwid <- c(0.1, 4.5)
margins <- c(0.5, 0.5)

# The minimum and maximum number of initial clusters to use
min.clusters <- 6
max.clusters <- 15

#min.clusters <- 9
#max.clusters <- 10

# Whether or not we should perform the 1D clustering along the margins.
doClusteringAlongMargins=FALSE
#doClusteringAlongMargins=TRUE

# Keep track of the total number of runs [as we vary number of clusters
# and the various methods (binomial, beta, and gaussian)]
total.num.runs <- 0

# C will the confusion/consensus clustering matrix summed across all
# methods (binomial, beta, and gaussian)
C <- 0

# C.binomial will the confusion/consensus clustering matrix for the binomial
# method
C.binomial <- 0

num.runs.per.method <- length(min.clusters:max.clusters)
col.palette <- bluered(num.runs.per.method)

# First exclude and points that we found in our previously clustering had low copy number,
# or bad copy number calls--these will have cluster NA.
# If we don't they will be consistently thrown out and will show up the consensus matrix as not being
# equal to themselves, since getConsensusMatrix sets entries to 0 if the cluster is NA.
v1 <- merge(v1, valid.results.chr.start, by.x=c(1,2), by.y=c(1,2), all.x=FALSE, all.y=FALSE)
v2 <- merge(v2, valid.results.chr.start, by.x=c(1,2), by.y=c(1,2), all.x=FALSE, all.y=FALSE)

# Do binomial clustering
cluster.with.binomial <- TRUE
if (cluster.with.binomial) {

  # Vary the number of initial clusters
  for(i in min.clusters:max.clusters) {  
    sc <- sciClone(vafs=list(v1,v2), sampleNames=samples, useSexChrs=FALSE, copyNumberCalls=list(cn1,cn2), copyNumberMargins=0.25, minimumDepth=100, clusterMethod="binomial.bmm", verbose=1, maximumClusters=i,doClusteringAlongMargins=doClusteringAlongMargins)

    # Write the VAF cluster assignments and the cluster summaries to a file    
    output <- paste("aml28.binomial.", i, sep="")
    writeClusterTable(sc, paste(output, ".clustered", sep=""))
    writeClusterSummaryTable(sc, paste(output, ".clusters", sep=""))

    # Output the 2D plot to a file
    plot.2d.file <- paste(output, ".2d.pdf", sep="")
    sc.plot2d(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, plot.2d.file)

    # Initial number clusters = min.clusters (i.e., 6) will cluster the
    # ambiguous points differently.
    # Save this file as supp1.
    if(i == 10) {
      system(paste("cp", plot.2d.file, "supp-fig2.pdf"))
    }
    
    # If we did 1D clustering along the margins, output that as well.    
    if(doClusteringAlongMargins) {
      sc.plot2dWithMargins(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, paste(output, ".2d.with-margins.pdf", sep=""))
    }

    # Get the connectivity/confusion/consensus matrix for this one run
    # and add it to the binomial consensus matrix.
    tmp <- getConnectivityMatrix(sc)
    C.binomial <- C.binomial + tmp

    # Write the matrix to a file    
    write.table(file=paste(output, ".C.tsv", sep=""), tmp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
    # Don't bother plotting this individual matrix
    #pdf(paste(output, ".C.pdf", sep=""))
    #image(tmp)
    #dev.off()

    total.num.runs <- total.num.runs + 1
  }

  # Write the summed matrix across all binomial runs to a file  
  output <- "aml28.binomial"
  write.table(file=paste(output, ".C.tsv", sep=""), C.binomial, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # Don't plot the matrix this way, it looks ugly.  
  #pdf(paste(output, ".C.pdf", sep=""))
  #image(C.binomial)
  #dev.off()

  # Plot the binomial consensus matrix as a heatmap and save to a file.
  heatmap.file <- paste(output, ".heatmap.pdf", sep="")
  pdf(heatmap.file)
  #my.heatmap.2(C.binomial, dendrogram="none", trace="none", labRow=rep("", dim(C.binomial)[1]), labCol=rep("", dim(C.binomial)[2]),key=FALSE,col=col.palette,lhei=lhei,lwid=lwid,margins=margins)
  #my.heatmap.2(C.binomial, dendrogram="none", trace="none", labRow=rep("", dim(C.binomial)[1]), labCol=rep("", dim(C.binomial)[2]),col=col.palette)
  my.heatmap.2(C.binomial / num.runs.per.method, dendrogram="none", trace="none", labRow=rep("", dim(C.binomial)[1]), labCol=rep("", dim(C.binomial)[2]), col=col.palette, lmat=lmat, lwid=lwid, lhei=lhei, density.info="none", margins=margins, revC=TRUE)
  dev.off()

}

# C.gaussian will the confusion/consensus clustering matrix for the gaussian
# method
C.gaussian <- 0

# Do gaussian clustering
cluster.with.gaussian <- TRUE
if (cluster.with.gaussian) {

  # Vary the number of initial clusters
  for(i in min.clusters:max.clusters) {
    sc <- sciClone(vafs=list(v1,v2), sampleNames=samples, useSexChrs=FALSE, copyNumberCalls=list(cn1,cn2), copyNumberMargins=0.25, minimumDepth=100, clusterMethod="gaussian.bmm", verbose=1, maximumClusters=i,doClusteringAlongMargins=doClusteringAlongMargins)

    # Write the VAF cluster assignments and the cluster summaries to a file
    output <- paste("aml28.gaussian.", i, sep="")
    writeClusterTable(sc, paste(output, ".clustered", sep=""))
    writeClusterSummaryTable(sc, paste(output, ".clusters", sep=""))

    plot.2d.file <- paste(output, ".2d.pdf", sep="")
    # Output the 2D plot to a file
    sc.plot2d(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, paste(output, ".2d.pdf", sep=""))

    if(i == 10) {
      system(paste("cp", plot.2d.file, "supp-fig1.pdf"))
    }

    
    # If we did 1D clustering along the margins, output that as well.
    if(doClusteringAlongMargins) {
      sc.plot2dWithMargins(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, paste(output, ".2d.with-margins.pdf", sep=""))
    }

    # Get the connectivity/confusion/consensus matrix for this one run
    # and add it to the gaussian consensus matrix.
    tmp <- getConnectivityMatrix(sc)
    C.gaussian <- C.gaussian + tmp

    # Write the matrix to a file
    write.table(file=paste(output, ".C.tsv", sep=""), tmp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Don't bother plotting this individual matrix
    #pdf(paste(output, ".C.pdf", sep=""))
    #image(tmp)
    #dev.off()

    total.num.runs <- total.num.runs + 1
  }

  # Write the summed matrix across all gaussian runs to a file
  output <- "aml28.gaussian"
  write.table(file=paste(output, ".C.tsv", sep=""), C.gaussian, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # Don't plot the matrix this way, it looks ugly.
  #pdf(paste(output, ".C.pdf", sep=""))
  #image(C.gaussian)
  #dev.off()

  # Plot the gaussian consensus matrix as a heatmap and save to a file.
  pdf(paste(output, ".heatmap.pdf", sep=""))
  #my.heatmap.2(C.gaussian, dendrogram="none", trace="none", labRow=rep("", dim(C.gaussian)[1]), labCol=rep("", dim(C.gaussian)[2]),key=FALSE,col=col.palette,lhei=lhei,lwid=lwid,margins=margins)
  my.heatmap.2(C.gaussian / num.runs.per.method, dendrogram="none", trace="none", labRow=rep("", dim(C.gaussian)[1]), labCol=rep("", dim(C.gaussian)[2]), col=col.palette, lmat=lmat, lwid=lwid, lhei=lhei, density.info="none", margins=margins, revC=TRUE)
  dev.off()
}

# C.beta will the confusion/consensus clustering matrix for the beta
# method
C.beta <- 0

# Do beta clustering
cluster.with.beta <- TRUE
if (cluster.with.beta) {

  # Vary the number of initial clusters  
  for(i in min.clusters:max.clusters) {  
    sc <- sciClone(vafs=list(v1,v2), sampleNames=samples, useSexChrs=FALSE, copyNumberCalls=list(cn1,cn2), copyNumberMargins=0.25, minimumDepth=100, verbose=1, maximumClusters=i,doClusteringAlongMargins=doClusteringAlongMargins)

    # Write the VAF cluster assignments and the cluster summaries to a file    
    output <- paste("aml28.beta.", i, sep="")
    writeClusterTable(sc, paste(output, ".clustered", sep=""))
    writeClusterSummaryTable(sc, paste(output, ".clusters", sep=""))

    # Output the 2D plot to a file    
    sc.plot2d(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, paste(output, ".2d.pdf", sep=""))

    # If we did 1D clustering along the margins, output that as well.    
    if(doClusteringAlongMargins) {
      sc.plot2dWithMargins(sc, highlightsHaveNames=highlightsHaveNames, positionsToHighlight=genesToAnnotate, paste(output, ".2d.with-margins.pdf", sep=""))
    }
    
    # Get the connectivity/confusion/consensus matrix for this one run
    # and add it to the beta consensus matrix.
    tmp <- getConnectivityMatrix(sc)
    C.beta <- C.beta + tmp

    # Write the matrix to a file    
    write.table(file=paste(output, ".C.tsv", sep=""), tmp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Don't bother plotting this individual matrix    
    #pdf(paste(output, ".C.pdf", sep=""))
    #image(tmp)
    #dev.off()

    total.num.runs <- total.num.runs + 1
  }

  # Write the summed matrix across all beta runs to a file  
  output <- "aml28.beta"
  write.table(file=paste(output, ".C.tsv", sep=""), C.beta, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # Don't plot the matrix this way, it looks ugly.  
  #pdf(paste(output, ".C.pdf", sep=""))
  #image(C.beta)
  #dev.off()

  # This is Fig 3A.
  # Plot the beta consensus matrix as a heatmap and save to a file.
  heatmap.file <- paste(output, ".heatmap.pdf", sep="")
  pdf(heatmap.file)
  #my.heatmap.2(C.beta, dendrogram="none", trace="none", labRow=rep("", dim(C.beta)[1]), labCol=rep("", dim(C.beta)[2]),key=FALSE,col=col.palette,lhei=lhei,lwid=lwid,margins=margins)
  my.heatmap.2(C.beta / num.runs.per.method, dendrogram="none", trace="none", labRow=rep("", dim(C.beta)[1]), labCol=rep("", dim(C.beta)[2]), col=col.palette, lmat=lmat, lwid=lwid, lhei=lhei, density.info="none", margins=margins, revC=TRUE)
  dev.off()
  #system(paste("cp", heatmap.file, "figure4a.pdf"))
}

# Make the "grand" consensus matrix C by summing the consensus matrices
# from all 3 methods.
C <- C.gaussian + C.binomial + C.beta

# Save this grand consensus matrix to a file
output <- "aml28"
write.table(file=paste(output, ".C.tsv", sep=""), C, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Ugly, don't plot like this.
#pdf(paste(output, ".C.pdf", sep=""))
#image(C)
#dev.off()

# Plot the grand consensus matrix as a heatmap and save to a file.
heatmap.file <- "aml28.heatmap.pdf"
pdf(heatmap.file)
col.palette <- bluered(total.num.runs)
#my.heatmap.2(C, dendrogram="none", trace="none", labRow=rep("", dim(C)[1]), labCol=rep("", dim(C)[2]),key=FALSE,col=col.palette,lhei=lhei,lwid=lwid,margins=margins)
my.heatmap.2(C / total.num.runs, dendrogram="none", trace="none", labRow=rep("", dim(C)[1]), labCol=rep("", dim(C)[2]), col=col.palette, lmat=lmat, lwid=lwid, lhei=lhei, density.info="none", margins=margins, revC=TRUE)
dev.off()
system(paste("cp", heatmap.file, "figure4a.pdf"))


# Plot Fig3B:  the number of clusters as we vary the number of VAFs

# Downsample by a fraction of the original number of points
orig.num.items <- dim(valid.results)[1]
ratio <- 0.75
num.ratios <- 16
num.iterations <- 10

# num.ratios <- 5
# num.iterations <- 2

# Also, store the number of clusters
nc <- matrix(data=0, nrow=num.ratios, ncol=num.iterations)

idx <- 1:orig.num.items

ratios <- rep(0, num.ratios)
down.sample.nums <- rep(0, num.ratios)
for(i in 1:num.ratios) {
  ratios[i] <- ratio^(num.ratios-i)
  down.sample.nums[i] <- round(orig.num.items * ratios[i])
  down.sample.num <- down.sample.nums[i]
  
  for(n in 1:num.iterations) {

    cat("ratio = ", ratios[i], "n = ", n, " num = ", down.sample.num, "\n")

    sample.indices <- sort(sample(idx, down.sample.num))
    
    # Randomly sample the original data
    random.subset <- valid.results.chr.start[sample.indices,]
    v1.new <- merge(v1, random.subset, by.x=c(1,2), by.y=c(1,2), all.x=FALSE, all.y=FALSE)
    v2.new <- merge(v2, random.subset, by.x=c(1,2), by.y=c(1,2), all.x=FALSE, all.y=FALSE)

    # Perform the clustering. NB: disable outlier detection (given that some
    # of these data sets will be very small).
    sc <- sciClone(vafs=list(v1.new,v2.new), sampleNames=samples, useSexChrs=FALSE, copyNumberCalls=list(cn1,cn2), copyNumberMargins=0.25, minimumDepth=100, verbose=1, doClusteringAlongMargins=FALSE, clusterParams="no.pvalue.outlier.detection")
    
    # Plot the 2D file
    sc.plot2d(sc, paste("fig4-r-", ratios[i], "-n-", n, ".pdf", sep=""))

    # Store the VAF cluster assignments.
    cluster.file <- paste("fig4-r-", ratios[i], "-n-", n, ".cluster", sep="")
    writeClusterTable(sc, cluster.file)

    # Read back in the VAF cluster assignments to find the number of clusters.
    # Should be a better way to do this.
    cluster.results.downsample <- read.table(cluster.file, sep="\t", header=TRUE)
    num.clusters <- length(unique(cluster.results.downsample$cluster))

    # Store the number of clusters
    nc[i,n] <- num.clusters
  }
  
}

# Calculate the average and std devs of the number of clusters
# and store them in a file.
nc.avg <- apply(nc, 1, mean)
nc.sd <- apply(nc, 1, sd)

tab <- data.frame(nc.avg=nc.avg, nc.sd=nc.sd, num.pts=down.sample.nums)
write.table(file="stats.tsv", tab, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Copy Chris' code instead, looks much better.
#pdf(file="figure3b.pdf", width=7.2, height=6, bg="white")
#plot(1:orig.num.items, ylim=c(0,num.clusters.orig), type="n", xlab="Num Variants", ylab="Num Clusters")
#points(down.sample.nums, nc.avg)
#plotCI(x=down.sample.nums, y=nc.avg, err="y", uiw=nc.sd, add=TRUE)
#dev.off()

pdf(file="figure4b.pdf", width=13, height=5, bg="white")
plot(1:850, ylim=c(0,6), type="n", xlab="Number of Downsampled Variants", ylab="Number of Clusters", xaxs="i", cex.axis=1.3, cex.lab=1.3)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey95")
abline(v=seq(0,850,50), col="grey50", lty=3)
abline(h=seq(0,6,1), col="grey50", lty=3)
points(down.sample.nums, nc.avg, type="b", col="darkgreen", pch=16, cex=1.3)
plotCI(x=down.sample.nums, y=nc.avg, err="y", uiw=nc.sd, add=TRUE, col="darkgreen", pch=16, cex=1.3)
dev.off()



save.image("out.Rdata")

#clean up the mess
dir.create("tmp")
files = c(Sys.glob("*.pdf"),Sys.glob("*.clust*"), Sys.glob("*.tsv"))
for(i in files){
  file.rename(i,paste("tmp",i,sep="/"))
}

#move back the files we want
file.rename("tmp/figure4a.pdf","figure4c.pdf")
file.rename("tmp/figure4b.pdf","figure4a.pdf")

file.rename("tmp/supp-fig1.pdf","suppFigure2a.pdf")
file.rename("tmp/supp-fig2.pdf","suppFigure2b.pdf")

