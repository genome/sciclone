####--------------------------------------------------------------------
## vafs is a data frame of somatic variants- chr,pos,refreads,varreads,vaf
##     vaf is a percentage
##     if you want to filter on minimum depth in normal or tumor, do that
##     before calling this function
##
## copyNumberCalls is a data frame - chr st sp cn
##
## positionsToHighlight is a data frame where the first two columns are - chr, pos
##
sciClone <- function(vafs, outputPrefix, copyNumberCalls=NULL, sampleNames,
                     minimumDepth=100, clusteredDataOutputFile=NULL,
                     clusterMethod="bmm",clusterParams=NULL,
                     minimumLabelledPeakHeight=0.001, onlyLabelHighestPeak=FALSE,
                     overlayClusters=FALSE, plotOnlyCN2=FALSE, positionsToHighlight=NULL,
                     purities=NULL, highlightSexChrs=TRUE, testing=FALSE, cnCallsAreLog2=FALSE,
                     useSexChrs=TRUE, highlightsHaveNames=FALSE, doClustering=TRUE, verbose=TRUE){



  if(verbose){print("checking input data...")}

  #how many samples do we have?
  dimensions=NULL;
  if(is.data.frame(vafs)){
    dimensions = 1;
    vafs=list(vafs)
    copyNumberCalls=list(copyNumberCalls)
    sampleNames=c(sampleNames)
  } else if(is.list(vafs)){
    dimensions = length(vafs);
  } else {
    stop("input param vafs must be either a data frame (for 1-sample clustering), or a list of data frames (for multi-sample clustering)")
  }

  ##some sanity checks on the input data
  if(dimensions > 1){
    if(length(sampleNames) != dimensions){
      stop(paste("the number of sample names (",length(sampleNames),") does not equal the number of input samples (",dimensions,")",sep=""))
    }
  }
  if(!(is.null(copyNumberCalls))){
    if(length(copyNumberCalls) != dimensions){
      stop(paste("the number of input copy number calls(",length(copyNumberCalls),") does not equal the number of input samples (",dimensions,")",sep=""))
    }
  } else {
    print("No copy number files specified. Assuming all variants have a CN of 2.")
    copyNumberCalls = vector("list",dimensions)
  }
  
  if(!(is.null(purities))){
    if(length(purities) != dimensions){
      stop(paste("the number of input purities calls(",length(purities),") does not equal the number of input samples (",dimensions,")\nEither provide a purity for each sample, or set purities to NULL and it will be estimated for you",sep=""))
    }
  }

  if(highlightsHaveNames){
    if(is.null(positionsToHighlight)){
      print("ERROR - if highlightsHaveNames is true, positionsToHighlight must be provided")
      stop()
xs    }
    plotOnlyCN2=TRUE;
  }
  
  
  
  ##-----------------------------------------
  if(verbose){print("calculating kernel density and purity")}

  densityData=NULL

  doPurityEst=FALSE;
  if(is.null(purities)){
    doPurityEst=TRUE;
    purities = c()
  }

  ##clean up data, get kernel density, estimate purity
  for(i in 1:dimensions){
    vafs[[i]] = cleanAndAddCN(vafs[[i]], copyNumberCalls[[i]], i, cnCallsAreLog2, useSexChrs, minimumDepth)
    ##calculate the densities and peaks for variants of each copy number
    if(is.null(densityData)){
      densityData = list(getDensity(vafs[[i]]))
    } else {
      densityData= c(densityData, list(getDensity(vafs[[i]])))
    }
    
    ## This stop should probably be replaced so that plotting can take place
    ## for tumors with ployplody, even if we can't cluster
    ## (maybe we should even cluster with 3x regions, etc - put it on the todo list)
    if(is.null(densityData[[i]]$densities[[2]])){
      stop(paste("can't do clustering - no copy number 2 regions to operate on in sample",i));
    }

    if(doPurityEst){
      purities = c(purities,getPurity(densityData[[i]]$peakPos))
    }
  }


  ##-----------------------------------------------
  ##do the clustering
  if(verbose){print("Doing clustering...")}

  ## merge the data frames to get a df with vafs and readcounts for each variant in each sample
  vafs.merged = vafs[[1]]

  if(dimensions > 1){
    vafs.merged = vafs[[i]];
    for(i in 2:dimensions){
      vafs.merged = merge(vafs.merged, vafs[[i]], by.x=c(1,2), by.y=c(1,2), suffixes=c(i-1,i), all.x=TRUE, all.y=TRUE)
    }
  }  
  
  refcols = grep("^ref",names(vafs.merged))
  varcols = grep("^var",names(vafs.merged))
  vafcols = grep("^vaf",names(vafs.merged))
  depthcols = grep("^depth",names(vafs.merged))
  cncols = grep("^cn",names(vafs.merged))
  ##change NA values introduced by merge to ref/var/vaf/depth of zero, cn of 2
  for(i in c(vafcols,refcols,varcols,depthcols)){
    vafs.merged[is.na(vafs.merged[,i]),i] = 0;
  }
  for(i in cncols){
    vafs.merged[is.na(vafs.merged[,i]),i] = 2;
  }
  
  #add sample names to make output pretty
  names(vafs.merged)[refcols] = paste(sampleNames,".ref",sep="")
  names(vafs.merged)[varcols] = paste(sampleNames,".var",sep="")
  names(vafs.merged)[vafcols] = paste(sampleNames,".vaf",sep="")
  names(vafs.merged)[depthcols] = paste(sampleNames,".depth",sep="")
  names(vafs.merged)[cncols] = paste(sampleNames,".cn",sep="")  
  
  ##---------------------------------------------------
  ##do the clustering
  
  ## remove any lines where all CN columns are not 2
  ## we only cluster based on sites that are CN neutral in all samples

  ## there is probably a better way to do this with an apply function...
  cnNeutral = rep(1,length(vafs.merged[,1]))
  for(i in 1:length(vafs.merged[,1])){
    if(sum(as.numeric(vafs.merged[i,cncols]==2)) < length(cncols)){
      cnNeutral[i] = 0
    }
  }
  vafs.merged.cn2 = vafs.merged[as.logical(cnNeutral),]

  #convert to a matrix to feed into clustering
  vafs.matrix = as.matrix(vafs.merged.cn2[,vafcols])
  #convert vafs to be between 0 and 1
  vafs.matrix = vafs.matrix/100

  clust=NULL
  if(doClustering){
    clust=clusterVafs(vafs.matrix, clusterMethod, purities, clusterParams)
  }
  numClusters=0
  if(!(is.null(clust))){
    numClusters = max(clust$cluster.assignments)
    #append cluster assignments
    vafs.merged.cn2 = cbind(vafs.merged.cn2,cluster=clust$cluster.assignments)
    vafs.merged = merge(vafs.merged,vafs.merged.cn2, by.x=c(1:length(vafs.merged)), by.y=c(1:length(vafs.merged)),all.x=TRUE)
    #sort by chr, st
    vafs.merged = vafs.merged[order(vafs.merged[,1], vafs.merged[,2]),]
    print(paste("found",numClusters,"clusters"))
  }

  ##--------------------------------------------------------
  ##output and plots
  if(!(is.null(clust)) & !(is.null(clusteredDataOutputFile))){
    ##TODO - order these first    
    write.table(vafs.merged, file=clusteredDataOutputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
  }

  ##--------------------------------------------------
  ##set up the plot
  plot1d(vafs.merged, outputPrefix, densityData, sampleNames, dimensions, plotOnlyCN2,
         clust, highlightSexChrs, positionsToHighlight, highlightsHaveNames,
         overlayClusters, onlyLabelHighestPeak, minimumLabelledPeakHeight);
#  plot2d(outputPrefix) #multiple times
}


##---------------------------------------------------------------------------------
## Create the one dimensional plot with kde and scatter
##
plot1d <- function(vafs.merged, outputPrefix, densityData, sampleNames, dimensions, plotOnlyCN2, clust, highlightSexChrs, positionsToHighlight, highlightsHaveNames, overlayClusters, onlyLabelHighestPeak, minimumLabelledPeakHeight){
  pdf(file=paste(outputPrefix,".1d.pdf",sep=""), width=3.3, height=7.5, bg="white");
  numClusters = max(clust$cluster.assignments)
  
  ##one plot for each sample
  for(d in 1:dimensions){
    par(mfcol=c(5,1),mar=c(0.5,3,1,1.5),oma=c(3,0,4,0),mgp = c(3,1,0));

    name=sampleNames[d]
    
    densities = densityData[[d]]$densities
    factors = densityData[[d]]$factors
    peakPos = densityData[[d]]$peakPos
    peakHeights = densityData[[d]]$peakHeights
    maxDepth = densityData[[d]]$maxDepth
    maxDensity = densityData[[d]]$maxDensity

  
    ##draw the density plot
    scalingFactor = 25/maxDensity;
    plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");

    ##plot bg color
    rect(0, 0, 100, 28, col = "#00000011",border=NA);
    
    axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxDensity)),las=1,
         cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);
    
    ##are we plotting everything or just CN2?
    cnToPlot = c();
    if(plotOnlyCN2){
      cnToPlot = c(2)
    } else {
      cnToPlot = 1:4
    }
        
    ##colors for different copy numbers
    colors=c("#1C3660AA","#67B32EAA","#F49819AA","#E52420AA")
    
    for(i in cnToPlot){
      if(!(is.null(densities[[i]]))){
        ##density lines
        lines(densities[[i]]$x, scalingFactor*factors[[i]], col=colors[i], lwd=2);
        ##peak labels
        ppos = c();
        if(onlyLabelHighestPeak){
          ppos = which(peakHeights[[i]] == max(peakHeights[[i]]))
        } else {
          ppos = which((peakHeights[[i]] == max(peakHeights[[i]])) &
            (peakHeights[[i]] > minimumLabelledPeakHeight))
        }
        if(length(ppos) > 0){
          text(x=peakPos[[i]][ppos], y=(scalingFactor * peakHeights[[i]][ppos])+1.7,
               labels=signif(peakPos[[i]][ppos],3),
               cex=0.7, srt=0, col=colors[[i]]);
        }
      }
    }
    
    
    ##if we have X/Y vals from the clustering alg, add them
    if(!(is.null(clust))){
      points(clust$fit.x, clust$fit.y*25, type="l",col="grey50")
    }
  
    ##legend
    leg = c("1 Copy","2 Copies","3 Copies","4 Copies")
    lcol=colors
    if( length(cnToPlot)== 1 ){
      leg=c("2 Copies")
      lcol=colors[2]
    }
    if(!(is.null(clust))){
      leg = c(leg,"Model Fit")
      lcol = c(lcol, "grey50")
    }
    legend(x="topright", lwd=2, legend=leg, col=lcol, bty="n", cex=0.6, y.intersp=1.25);
    
    
    axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
    mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
    mtext("Kernel Density",side=2,cex=0.5,padj=-4.2);
    
    ##add a title to the plot
    title=""
    if(is.null(sampleNames)){
      title="Clonality Plot"
    } else {
      title=paste(sampleNames,"Clonality Plot",sep=" ");
    }
    mtext(title, adj=0.5, padj=-5, cex=0.65, side=3);
    

    ##-----------------------------------------------------
    ##create the scatterplots of vaf vs density

    ##grab only the vafs for this sample:
    a = vafs.merged[,c("chr","st")]
    n = (d*5)-2 #this sample's columns
    b = vafs.merged[ ,n:(n+4)]
    vafs = cbind(a,b)
    names(vafs) = c("chr","st","ref","var","vaf","depth","cn")
    if(numClusters>0){
      vafs = cbind(vafs, vafs.merged$cluster)
      names(vafs)[8] = "cluster"
    }


    for(i in cnToPlot){
      v = vafs[vafs$cn==i,];
      print(length(v[,1]))
      drawScatterPlot(v, highlightSexChrs, positionsToHighlight, colors, i, maxDepth, highlightsHaveNames, overlayClusters)
      axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.4);
    }
  }
    
  
  if(length(cnToPlot < 2) && highlightsHaveNames){
    addHighlightLegend(vafsByCn[[i]], positionsToHighlight)
  } else {
    if(highlightsHaveNames){
      print("WARNING: highlighted point naming is only supported when plotOnlyCN2 is TRUE")
    }
  }
  
  
  ##close the pdf
  devoff <- dev.off();
}


##--------------------------------------------------------------------
## intersect the variants with CN calls to classify them
##
addCnToVafs <- function(vafs,cncalls){
  library(IRanges)  
  vafs$cn = NA
  ##for each chromosome
  for(chr in names(table(vafs$chr))){
    vars = IRanges(start=vafs[vafs$chr==chr,]$st, end=vafs[vafs$chr==chr,]$st)
    cnsegs = IRanges(start=cncalls[cncalls[,1]==chr,2], end=cncalls[cncalls[,1]==chr,3])
    if( (length(vars) == 0) | (length(cnsegs) == 0)){
      next
    }
    matches = as.matrix(findOverlaps(vars,cnsegs))
    if(length(matches) > 0){
      for(i in 1:length(vars)){
        vpos = which(vafs$chr==chr & vafs$st==start(vars[matches[which(matches[,1]==i),1]]))
        if(identical(vpos,integer(0))) { next }
        cpos = which(cncalls[,1] == chr & cncalls[,2]==start(cnsegs[matches[which(matches[,1]==i),2]]))
        ##set the value
        vafs[vpos,]$cn = cncalls[cpos,4]
      }
    }
  }

  if(length(which(is.na(vafs$cn))) > 0){
    print("Not all variants fall within a provided copy number region. The copy number of these variants is assumed to be 2.")
    vafs[which(is.na(vafs$cn)),]$cn = 2
  }
  
  return(vafs)
}



##---------------------------------------------------------------------
## clean up vaf data, add cn
##
cleanAndAddCN <- function(vafs, cn, num, cnCallsAreLog2, useSexChrs, minimumDepth){
  names(vafs) = c("chr","st","ref","var","vaf")

  ##remove MT values
  vafs = vafs[!(vafs$chr == "M" | vafs$chr == "MT"),]

  ##remove NA sites
  vafs = vafs[!(is.na(vafs$vaf)),]

  ##add depth
  vafs = vafs[vafs$vaf > 0,]
  vafs$depth = round(vafs$var/(vafs$vaf/100))

  ##add cn calls
  if(is.null(cn)){
    ##assume all sites are 2x if no cn info
    vafs$cn = 2;
  } else {
    if(cnCallsAreLog2){
      cn[,4] = (2^(copyNumberCalls[,4]))*2
    }
    vafs = addCnToVafs(vafs,cn)
  }
  ##remove sex chromosomes if specified
  if(!(useSexChrs)){
    vafs = vafs[vafs$chr != "X" & vafs$chr != "Y",]
  }
  
  ##remove any sites with less than the minimum depth
  vafs = vafs[vafs$depth >= minimumDepth,]
  if(length(vafs$chr) == 0){
    print(paste("No variants in sample",num,"exceed a depth of",minimumDepth,". Lower this threshold and try again."))
    stop()
  }
  print(paste("Number of variants with depth >= ",minimumDepth," in sample ",num," being used for analysis: ",length(vafs$chr),sep=""))

  return(vafs)

}



##--------------------------------------------------------------------------
## calculate a samples purity from VAF peaks
##
getPurity <- function(peakPos){
  purity = 0
  if(length(peakPos[[2]][peakPos[[2]] <= 60]) > 0){
    purity = max(peakPos[[2]][peakPos[[2]] <= 50])*2;
    if(length(min(peakPos[[2]][peakPos[[2]] > 50])) > 0){
      nextHighestPeak = min(peakPos[[2]][peakPos[[2]] > 50]);
      ##if the peak is between 50 and 60, assume that it's noise and
      ##the peak is actually at 50
      if (nextHighestPeak > 50 && nextHighestPeak < 60 && purity < 60) {
        purity = 100;
      }
    }
  }
  ##if all of these failed to find a good peak, assume 100% purity
  if (purity == 0) {
    print("Unable to make reliable calculation of purity, so assuming 100%")
    purity = 100;
  }
  print(paste("Tumor purity estimated to be ",signif(purity,4),".",sep=""));
  return(purity)
}



##--------------------------------------------------------------------------
## calculate the 1d kernel density and peaks
##
getDensity <- function(vafs){
  ##data structure to hold infoemacs sp
  densities = vector("list",4)
  factors = vector("list",4)
  peakPos = vector("list",4)
  peakHeights = vector("list",4)
  maxDensity = 0
  maxDepth=0

  ##default cutoffs are +/- 0.5x cn
  cnLoThresh = c(0,1.5,2.5,3.5)
  cnHiThresh = c(1.5,2.5,3.5,4.5)

  for(i in 1:4){
    ##grab only the variants in this copy number
    v = vafs[vafs$cn > cnLoThresh[i] & vafs$cn < cnHiThresh[i],]
    if(length(v[,1]) > 0){

      ##need two points for density calc
      if(length(v[,1]) > 1){
        ##calculate the density
        densities[[i]] = density(v$vaf, from=0, to=100, na.rm=TRUE)
        factors[[i]] = (length(v[,1])/length(vafs[,1]))*densities[[i]]$y
        ##find the peaks
        p = c(getPeaks(factors[[i]]),FALSE,FALSE)
        peakPos[[i]] = densities[[i]]$x[p]
        peakHeights[[i]] = factors[[i]][p]
        ##store the largest density for use in scaling the plots later
        if(max(factors[[i]]) > maxDensity){
          maxDensity = max(factors[[i]])
        }

      }

      ##store the largest depth for use in scaling the plots later
      if(max(v$depth) > maxDepth){
        maxDepth = max(v$depth)
      }

    } #else has a value of NULL
  }
  return(list(densities=densities,
              factors=factors,
              peakPos=peakPos,
              peakHeights=peakHeights,
              maxDensity=maxDensity,
              maxDepth=maxDepth))
}





##--------------------------------------------------------------------
## find inflection points (peaks)
##
getPeaks<-function(series,span=3){
  z <- embed(series, span);
  s <- span%/%2;
  v<- max.col(z) == 1 + s;
  result <- c(rep(FALSE,s),v);
  result <- result[1:(length(result)-s)];
  return(result)
}



##--------------------------------------------------------------------
## draw a scatter plot of vaf vs depth
##
drawScatterPlot <- function(data, highlightSexChrs, positionsToHighlight, colors, cn, maxDepth, highlightsHaveNames, overlayClusters){

  ## define plot colors
  ptcolor = colors[cn]
  circlecolor = substr(colors[cn],1,7) #chop off the alpha value

  ## define the plot space by plotting offscreen points
  plot.default( x=-10000, y=1, log="y", type="p", pch=19, cex=0.4,
               col="#00000000", xlim=c(0,100), ylim=c(5,maxDepth*3),
               axes=FALSE, ann=FALSE, xaxs="i", yaxs="i");

  addPoints <- function(data, color, highlightSexChrs, pch=NULL){
    outlineCol = rgb(0,0,0,0.1);
    if(highlightSexChrs){
      ##plot autosomes
      data.autosomes = data[!(data$chr == "X" | data$chr == "Y"),]
      points(data.autosomes$vaf, data.autosomes$depth, type="p", pch=16, cex=0.75, col=color);
      #points(data.autosomes$vaf, data.autosomes$depth, type="p", pch=1, cex=0.75, col=outlineCol, lwd=);
      ##plot sex chromsomes with different shape
      data.sex = data[(data$chr == "X" | data$chr == "Y"),]
      points(data.sex$vaf, data.sex$depth, type="p", pch=17, cex=0.75, col=color);
      #points(data.sex$vaf, data.sex$depth, type="p", pch=1, cex=0.75, col=outlineCol);
    } else {
      points(data$vaf, data$depth, type="p", pch=16, cex=0.75, col=color);
      ##add outline
      #points(data$vaf, data$depth, type="p", pch=1, cex=0.75, col=outlineCol);
    }

  }

  #do we have any points to plot?
  if(length(data[,1]) > 0){
    ##if we have cluster assignments in col 8, color them
    if(length(data) > 7 & overlayClusters){
      library(RColorBrewer)
      clusters=sort(unique(data[,8]))
      cols = NA;
      ##repeat colors if more than 5 clusters
      ## cols=brewer.pal(5,"Greens")
      ## if(length(clusters)>5){
      ##   cols = rep(ceiling(length(clusters)/5),cols)
      ## }
      cols=brewer.pal(8,"Dark2")
      if(length(clusters)>8){
        cols = rep(ceiling(length(clusters)/8),cols)
      }

      #add transparency to colors
      for(i in 1:length(cols)){
        z = col2rgb(cols[i])
        cols[i] = rgb(z[1], z[2], z[3], 200, maxColorValue=255)
      }

      for(i in 1:length(clusters)){
        p = data[data$cluster == i,]
        addPoints(p,cols[i],highlightSexChrs)
      }
    } else { #just use the normal color
      addPoints(data,ptcolor,highlightSexChrs)
    }

    ##TODO - add a legend for point types - highlight vs sex chrs vs autosomes

    ## add highlighted points selected for by user
    if(!(is.null(positionsToHighlight))){
      addpts = merge(data, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))
      if(length(addpts[,1]) > 1){
        if(highlightsHaveNames){
          for(i in 1:length(addpts$vaf)){
            text(addpts$vaf[i],addpts$depth[i],labels=i,cex=0.5)
          }
        } else {
          addPoints(addpts, col="#555555FF", highlightSexChrs);
          #points(x=addpts$vaf,y=addpts$depth,type="p",pch=7,cex=0.8,col="#555555FF");
        }
      }
    }
  }

  ## define the axis
  axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
  for (i in 2:length(axTicks(2)-1)) {
    lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
  }

  ## plot the background color
  rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA);

  ## add cn circle
  points(x=c(97),y=c(maxDepth),type="p",pch=19,cex=3,col=circlecolor);
  text(c(97),y=c(maxDepth), labels=c(cn), cex=1, col="#FFFFFFFF")


  ## y axis label
  mtext("Tumor Coverage",side=2,cex=0.5,padj=-4.2);
}


##--------------------------------------------------------------------
## add a legend for highlighted points with names
addHighlightLegend <- function(data, positionsToHighlight){
  plot.default(x=-10000, y=1, type="p", pch=19, cex=0.4,
               col="#00000000", xlim=c(0,1000), ylim=c(0,1000),
               axes=FALSE, ann=FALSE);

  mtext("Tier 1 Genes",side=2,cex=0.5,padj=-4.2);


  names(positionsToHighlight)=c("chr","st","name");

  if(!(is.null(positionsToHighlight))){
    addpts = merge(data, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))
    if(length(addpts[,1]) > 1){
      ypos=rev(seq(0,900,(900/13)))[1:13]
      ncol=ceiling(length(addpts[,1])/13)
      xpos=0;
      offset=1
      for(i in 1:ncol){
        names = addpts[offset:(offset+12),]$name;
        names = as.character(names[!(is.na(names))])
        num = length(names)

        for(i in 1:num){
          text(xpos, ypos[i], paste(offset+i-1,". ",names[i],sep=""), cex=0.5, pos=4)
        }
        xpos=xpos+250;
        offset=offset+13;
      }
    }
  }
}
