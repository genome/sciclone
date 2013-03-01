
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
      points(clust$fit.x, clust$fit.y[d,]*25, type="l",col="grey50")
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
    if(is.null(sampleNames[d])){
      title="Clonality Plot"
    } else {
      title=paste(sampleNames[d],"Clonality Plot",sep=" ");
    }
    mtext(title, adj=0.5, padj=-5, cex=0.65, side=3);


    ##-----------------------------------------------------
    ##create the scatterplots of vaf vs density

    ##grab only the vafs for this sample:
    vafs = getOneSampleVafs(vafs.merged, d, numClusters);

    for(i in cnToPlot){
      v = vafs[vafs$cn==i,];
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
      numClusters=length(unique(data[,8]))
      cols = getClusterColors(numClusters)
      for(i in 1:numClusters){
        p = data[data$cluster == i,]
        addPoints(p,cols[i],highlightSexChrs,pch=i)
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
          ##points(x=addpts$vaf,y=addpts$depth,type="p",pch=7,cex=0.8,col="#555555FF");
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


##---------------------------------------------------------------------------------
## Create the one dimensional plot with kde and scatter
##
plot2d <- function(vafs.merged, outputPrefix, sampleNames, dimensions, positionsToHighlight, highlightsHaveNames, overlayClusters){
  pdf(file=paste(outputPrefix,".2d.pdf",sep=""), width=7.2, height=6, bg="white")
  numClusters = max(vafs.merged$cluster)

  
  ##create a 2d plot for each pairwise combination of samples
  for(d1 in 1:(dimensions-1)){
    for(d2 in d1:dimensions){
      if(d1==d2){
        next
      }
      
      vafs1 = getOneSampleVafs(vafs.merged, d1, numClusters);
      vafs2 = getOneSampleVafs(vafs.merged, d2, numClusters);

      ##get only cn2 points
      vafs1 = vafs1[vafs1$cn==2,]
      vafs2 = vafs2[vafs2$cn==2,]

      v = merge(vafs1,vafs2,by.x=c(1,2,8), by.y=c(1,2,8),suffixes=c(".1",".2"))

      cols = getClusterColors(numClusters)
      #create the plot
      #layout(matrix(c(1,2),1,2, byrow=TRUE), widths=c(4,1), heights=c(1,1))
      plot(-100, -100, xlim=c(0,120), ylim=c(0,100), main=paste(sampleNames[d1],"vs",sampleNames[d2]),
           xlab=paste(sampleNames[d1],"VAF                   "), ylab=paste(sampleNames[d2],"VAF"),
           bty="n", xaxt="n")
      abline(v=seq(0,100,20),col="grey50", lty=3)

      segments(rep(-10,5),seq(0,100,20),rep(105,5),seq(0,100,20), lty=3, col="grey50")
               
      #abline(h=seq(0,100,20),col="grey50", lty=3)

      axis(side=1,at=seq(0,100,20),labels=seq(0,100,20))

      
      for(i in 1:numClusters){
        if(overlayClusters){
          points(v[v$cluster==i,]$vaf.1, v[v$cluster==i,]$vaf.2, col=cols[i], pch=i)
        } else {
          points(v[v$cluster==i,]$vaf.1, v[v$cluster==i,]$vaf.2, pch=14)
        }
      }
      #plot(-100, -100, xlim=c(0,100), ylim=c(0,100),main="Clusters")
      legend("topright", legend=1:numClusters, col=cols[1:numClusters], title="Clusters", pch=1:numClusters)
    }
  }
  devoff = dev.off()
}



##-------------------------------------------------------------------------------------
## Get a list of colors to use for the clusters (after 20 they start repeating)
##
getClusterColors <- function(numClusters){
  library(RColorBrewer)
  cols = NA;

  cols=brewer.pal(8,"Dark2")
  if(numClusters>8){
    cols = c(cols,brewer.pal(9,"Set1"))
  }
  if(numClusters > 20){ #if this happens, something is probably weird with your data, but we'll support it anyway
    cols = rep(cols,ceiling(length(cols)/20))
  }
  cols = cols[1:numClusters]
  
  ##add transparency to colors
  for(i in 1:length(cols)){
    z = col2rgb(cols[i])
    cols[i] = rgb(z[1], z[2], z[3], 200, maxColorValue=255)
  }
  return(cols)  
}


##-------------------------------------------------------------------------------------
## Extract a single sample's data from the merged data frame
##
getOneSampleVafs <- function(vafs.merged, d, numClusters){
  a = vafs.merged[,c("chr","st")]
  n = (d*5)-2 #this sample's columns
  b = vafs.merged[ ,n:(n+4)]
  vafs = cbind(a,b)
  names(vafs) = c("chr","st","ref","var","vaf","depth","cn")
  if(numClusters>0){
    vafs = cbind(vafs, vafs.merged$cluster)
    names(vafs)[8] = "cluster"
  }
  return(vafs)
}
  
