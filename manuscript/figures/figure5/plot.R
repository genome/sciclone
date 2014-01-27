#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(sciClone))

load("out.Rdata")
sc.plot2d(sc,"figure5abc.pdf", singlePage=TRUE, scale=1.8)

## #do this by hand - 3d plotting doesn't work without interactive session
## sco = sc
## samplesToPlot=sc@sampleNames
## size=700
## open3d(windowRect = c(0,0, size, size) )

## a = sco@vafs.merged[,c(paste(samplesToPlot,".vaf",sep=""),"cluster")]
## a = a[!is.na(a$cluster),]
## a = a[!(a$cluster==0),]
## numClusters=max(a$cluster)
## cols=getClusterColors(numClusters)
## colvec = cols[a$cluster]
## plot3d(a[,1], a[,2], a[,3], xlim=c(0,100), ylim=c(0,100),zlim=c(0,100), axes=FALSE,
##        xlab="",ylab="", zlab="",
##        type="s", col=colvec)
## axes3d( edges=c("x--", "y--", "z"),labels=FALSE)
## for(i in c("+","-")){
##   for(j in c("+","-")){
##     axes3d( edges=paste("x",i,j,sep=""), tick=FALSE, labels=FALSE)
##     axes3d( edges=paste("y",i,j,sep=""), tick=FALSE, labels=FALSE)
##     axes3d( edges=paste("z",i,j,sep=""), tick=FALSE, labels=FALSE)
##   }
## }
## #position the plot the way you like it first, then do this to screencap
## rgl.postscript( "figure5d.pdf", fmt="pdf")
