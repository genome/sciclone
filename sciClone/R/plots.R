#---------------------------------------------------------------------------------
## Create the one dimensional plot with kde and scatter
##
sc.plot1d <- function(sco, outputFile,
                   plotOnlyCN2=FALSE, showCopyNumberScatterPlots=TRUE, highlightSexChrs=TRUE,
                   positionsToHighlight=NULL, highlightsHaveNames=FALSE, overlayClusters=TRUE,
                   overlayIndividualModels=TRUE, show1DHistogram=FALSE, onlyLabelHighestPeak=FALSE,
                   minimumLabelledPeakHeight=0.001, showTitle=TRUE){


  #sanity checks
  if(highlightsHaveNames){
    if(!(is.null(positionsToHighlight))){
      if(length(positionsToHighlight) < 3){
        print("ERROR: named plot requires names in the third column of the positionsToHighlight data frame")
        return(0)
      }
      plotOnlyCN2=TRUE
    } else {
      print("ERROR: highlightsHaveNames requires a 3-column dataframe of positions and names (chr, pos, name)");
      return(0);
    }
  }

  densityData = sco@densities
  vafs.merged = sco@vafs.merged

  sampleNames = sco@sampleNames
  dimensions = sco@dimensions

  clust = NULL
  if(overlayClusters){
    if(is.null(sco@clust[1])){
      print("ERROR: can't overlay clusters when clustering was not done on the input data")
      return(0)
    } else {
      clust = sco@clust
    }
  }

  pdf(file=outputFile, width=3.3, height=7.5, bg="white");

  numClusters = 0
  if(!is.null(clust)) {
    numClusters = max(clust$cluster.assignments)
  }

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
    scalingFactor = 1/maxDensity;
    plot.default(x=c(1:10),y=c(1:10),ylim=c(0,1.1), xlim=c(0,100), axes=FALSE,
                 ann=FALSE, col="#00000000", xaxs="i", yaxs="i");

    ##plot bg color
    rect(0, 0, 100, 1.1, col = "#00000011",border=NA);
    axis(side=2,at=c(0,1),labels=c(0, ""), las=1, cex.axis=0.6, hadj=0.6,
         lwd=0.5, lwd.ticks=0.5, tck=-0.01);


    ##are we plotting everything or just CN2?
    cnToPlot = c();
    if(plotOnlyCN2){
      cnToPlot = c(2)
    } else {
      cnToPlot = 1:4
    }
    ##grab only the vafs for this sample:
    vafs = getOneSampleVafs(vafs.merged, d, numClusters);
    ##colors for different copy numbers
    colors=c("#1C3660AA","#67B32EAA","#F49819AA","#E52420AA")

    for(i in cnToPlot){
      if(!(is.null(densities[[i]])) & (!show1DHistogram | (i!= 2))){
        ##density lines
        lines(densities[[i]]$x, scalingFactor*factors[[i]], col=colors[i], lwd=2);
        ##peak labels
        if(length(peakHeights[[i]]) > 0){
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
      } else if(show1DHistogram & (i == 2)) {

        ## Only show histogram for copy number = 2
        v = vafs[vafs$cleancn==2,];

        frequencies <- data.frame(x=v$vaf, row.names=NULL, stringsAsFactors=NULL)
        bin.width <- 2.5
        num.breaks <- ceiling(100/bin.width) + 1
        breaks <- unlist(lapply(0:(num.breaks-1), function(x) 100*x/(num.breaks-1)))
        h <- hist(v$vaf, breaks=breaks, plot=FALSE)
        ## Rescale the intensity so it has max value y = 1.
        h$density <- h$density / max(h$density)
        plot(h, add=TRUE, freq=FALSE, col="white", border="black")
      }
    }



    ##if we have X/Y vals from the clustering alg, add them
    ## If we overlay the model, show it in a different style
    ## style 4 = dot dashed
    ## style 1 = line
    model.style <- 4
    model.style <- 1
    model.width <- 0.25
    ## Plot the individual models with dotted lines (3) or dashed (2)
    individual.model.style <- 3
    individual.model.width <- 1

    if(!(is.null(clust))){
      maxFitDensity <- max(clust$fit.y[d,])
      #points(clust$fit.x, clust$fit.y[d,]*25, type="l",col="grey50")
      lines(clust$fit.x, clust$fit.y[d,]/maxFitDensity, type="l",col="grey50",lty=model.style, lwd=model.width)
      if(overlayIndividualModels==TRUE) {
        for(i in 1:numClusters) {
          lines(clust$fit.x, clust$individual.fits.y[[i]][d,]/maxFitDensity,
                type="l",col="grey50",lty=individual.model.style, lwd=individual.model.width)
        }
      }
    }

    ##legend
    lcol=colors
    lty = c(1,1,1,1)
    lwd = c(2,2,2,2)
    pchs = c(NA, NA, NA, NA)
    pt.bgs = lcol
    leg = c("1 Copy", "2 Copies", "3 Copies", "4 Copies")
    if( (length(cnToPlot)== 1) ) {
      if(show1DHistogram == FALSE) {
        lcol=colors[2]
        lty = c(1)
        lwd = c(2)
        pchs = c(NA)
        pt.bgs = lcol
        leg = c("2 Copies")
      } else {
        lcol="black"
        lty = c(0)
        lwd = c(0)
        pchs = c(22)
        pt.bgs = "white"
        leg = c("2 Copies")
      }
    }
    if(!(is.null(clust))){
      leg = c(leg,"Model Fit")
      lcol = c(lcol, "grey50")
      lty = c(lty, model.style)
      lwd = c(lwd, model.width)
      pt.bgs = c(pt.bgs, "grey50")
      pchs = c(pchs, NA)
      if(overlayIndividualModels==TRUE) {
        leg = c(leg,"Component Fits")
        lcol = c(lcol, "grey50")
        lty = c(lty, individual.model.style)
        lwd = c(lwd, 2)
        pt.bgs = c(pt.bgs, "grey50")
        pchs = c(pchs, NA)
      }
    }
    legend(x="topright", lwd=lwd, lty=lty, legend=leg, col=lcol, bty="n", cex=0.6, y.intersp=1.25, pch=pchs, pt.bg = pt.bgs);


    axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
    mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
    mtext("Density (a.u.)",side=2,cex=0.5,padj=-4.2);


    ##add a title to the plot
    if(showTitle){
      title=""
      if(is.null(sampleNames[d])){
        title="Clonality Plot"
      } else {
        title=paste(sampleNames[d],"Clonality Plot",sep=" ");
      }
      mtext(title, adj=0.5, padj=-5, cex=0.65, side=3);
    }

    ##-----------------------------------------------------
    ##create the scatterplots of vaf vs copy number

    if(showCopyNumberScatterPlots) {
      for(i in cnToPlot){
        v = vafs[vafs$cleancn==i,];
        drawScatterPlot(v, highlightSexChrs, positionsToHighlight, colors, i, maxDepth, highlightsHaveNames, overlayClusters)
        axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.4);


        if(length(cnToPlot) < 2 & highlightsHaveNames){
          addHighlightLegend(v, positionsToHighlight)
        } else {
          if(highlightsHaveNames){
            print("WARNING: highlighted point naming is only supported when plotOnlyCN2 is TRUE")
          }
        }
      }
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
    if(length(data) > 7 & overlayClusters & cn==2){
      numClusters=max(data[,8],na.rm=T)
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


##---------------------------------------------------------------------
## create the two dimensional plot with scatter annotated with
## clustering results and 1D plots along margins, this time using
## ggplot2

plot2dWithMargins <- function(sco, outputFile,positionsToHighlight=NULL, highlightsHaveNames=FALSE) {
  densityData = sco@densities
  vafs.merged = sco@vafs.merged
  vafs.1d.merged = sco@vafs.1d
  sampleNames = sco@sampleNames
  dimensions = sco@dimensions
  clust = sco@clust
  marginalClust = sco@marginalClust

  library(grid)
  library(ggplot2)

  plots.1d.list <- list()
  res.1d.max.densities <- list()

  xmin <- -5
  xmax <- 105

  tmp.file <- tempfile("outputFile.tmp")
  pdf(file=tmp.file, width=7.2, height=6, bg="white")

  # Create (and store) all possible 1D plots
  for(d in 1:dimensions){
    # Overlay the histogram of the data on the model fit.
    #ylab <- "\nDensity\n"
    ylab <- "Density (a.u.)"

    title <- ""

    # Set max posterior density to max of splinefun
    limits <- data.frame(x=c(min(marginalClust[[d]]$fit.x), max(marginalClust[[d]]$fit.x)))
    f <- splinefun(marginalClust[[d]]$fit.x, marginalClust[[d]]$fit.y[1,,drop=FALSE])
    max.posterior.density <- max(unlist(lapply(seq(from=limits$x[1], to=limits$x[2], by=10^-3), f)))

    # xlab <- paste("\n", sampleNames[d], "\n", sep="")
    xlab <- paste(sampleNames[d], " VAF", sep="")

    numClusters = 0
    if(!is.null(vafs.1d.merged[[d]]$cluster)) {
      numClusters = max(vafs.1d.merged[[d]]$cluster, na.rm=T)
    }

    vafs = getOneSampleVafs(vafs.1d.merged[[d]], d, numClusters)

    # Only show copy number = 2
    v = vafs[vafs$cleancn==2,];

    frequencies <- data.frame(x=v$vaf, row.names=NULL, stringsAsFactors=NULL)

    g <- ggplot(data = frequencies, aes(x)) + ggtitle(title) + xlab(xlab) + ylab(ylab)

    # g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

    g <- g + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

    g <- g + theme(plot.margin = unit(c(0,0,0,0), "cm"))

    bin.width <- 2.5
    num.breaks <- ceiling(100/bin.width) + 1
    breaks <- unlist(lapply(0:(num.breaks-1), function(x) 100*x/(num.breaks-1)))

    g <- g + geom_histogram(data = frequencies, mapping=aes(x, y=..ncount..*100), fill="white", colour="black", breaks=breaks)

    # Need to "print" the graph in order to see its maximum y value
    # NB:  this is redundant at this point, given that I scale the
    # histogram to have height 100 above.  i.e., max.density will be 100.
    tmp <- print(g)
    max.density <- max(tmp[["data"]][[1]]$ymax)
    res.1d.max.densities[[d]] <- max.density

    hline <- data.frame(x = c(0,100), y=c(-5,-5))
    g <- g + geom_line(data = hline, aes(x,y))

    #vline <- data.frame(y = c(0,max.density * 1.1), x=c(-5,-5))
    vline <- data.frame(y = c(0,100), x=c(-5,-5))
    g <- g + geom_line(data = vline, aes(x,y))

    scale <- max.density / max.posterior.density

    f <- splinefun(marginalClust[[d]]$fit.x, scale*marginalClust[[d]]$fit.y[1,,drop=FALSE])

    g <- g + stat_function(data = limits, fun=f, mapping=aes(x))

    g <- g + coord_cartesian(ylim=c(xmin, max.density*1.1), xlim=c(xmin,xmax))

    plots.1d.list[[d]] <- g
  }

  devoff = dev.off()
  unlink(tmp.file)
  
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip="off")

  pdf(file=outputFile, width=7.2, height=6, bg="white")

  ##create a 2d plot for each pairwise combination of samples
  for(d1 in 1:(dimensions-1)){
    for(d2 in d1:dimensions){
      if(d1==d2){
        next
      }

      xlab <- sampleNames[d1]
      xlab <- gsub("\\.", " ", xlab)
      #xlab <- paste("\n", xlab, "\n", sep="")
      xlab <- paste(xlab, " VAF", sep="")

      ylab <- sampleNames[d2]
      ylab <- gsub("\\.", " ", ylab)
      #ylab <- paste("\n", ylab, "\n", sep="")
      ylab <- paste(ylab, " VAF", sep="")

      numClusters = 0
      if(!is.null(vafs.merged$cluster)) {
        numClusters = max(vafs.merged$cluster, na.rm=T)
      }

      vafs1 = getOneSampleVafs(vafs.merged, d1, numClusters);
      vafs2 = getOneSampleVafs(vafs.merged, d2, numClusters);

      ##get only cn2 points
      vafs1 = vafs1[vafs1$cleancn==2,]
      vafs2 = vafs2[vafs2$cleancn==2,]

      v = merge(vafs1,vafs2,by.x=c(1,2,8), by.y=c(1,2,8),suffixes=c(".1",".2"))

      v.no.highlight <- v
      if(!(is.null(positionsToHighlight))) {
        chr.start.v <- cbind(v[,1], v[,2])
        chr.start.highlight <- cbind(positionsToHighlight[,1], positionsToHighlight[,2])        
        v.no.highlight <- v[!(apply(chr.start.v, 1, paste, collapse="$$") %in% apply(chr.start.highlight, 1, paste, collapse="$$")),]

      }

      title <- ""

      clusters <- v.no.highlight$cluster
      
      # Plot the points that we will not highlight
      frequencies.no.highlight <- data.frame(x=v.no.highlight$vaf.1, y=v.no.highlight$vaf.2, row.names=NULL, stringsAsFactors=NULL)

      g <- ggplot(data = frequencies.no.highlight, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = frequencies.no.highlight, aes(x=x, y=y), shape=clusters, colour=clusters) 

      # Now overlay any points that we will highlight
      if(!(is.null(positionsToHighlight))) {
        # Merge the data and the positions to highlight by chr (col 1)
        # and start (col 2)
        addpts = merge(v, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))        
        frequencies.highlight <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, row.names=NULL, stringsAsFactors=NULL)

        g <- g + geom_point(data = frequencies.highlight, aes(x=x, y=y), shape="*", size=10, colour="black")

      }
      
      # g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

      g <- g + theme_bw() + theme(panel.border = element_blank())

      g <- g + theme(plot.margin = unit(c(0,0,0,0), "cm"))

      hline <- data.frame(x = c(0,100), y=c(-5,-5))
      g <- g + geom_line(data = hline, aes(x,y))

      vline <- data.frame(y = c(0,100), x=c(-5,-5))
      g <- g + geom_line(data = vline, aes(x,y))

      g <- g + coord_cartesian(xlim=c(xmin, xmax), ylim=c(xmin, xmax))
      
      plot.2d <- g

      # Code to override clipping
      #plot.2d <- ggplot_gtable(ggplot_build(g))
      #plot.2d$layout$clip[plot.2d$layout$name=="panel"] <- "off"
      #grid.draw(plot.2d)

      plot.1d.1 <- plots.1d.list[[d1]]
      plot.1d.2 <- plots.1d.list[[d2]]

      grid.newpage()

      pushViewport(viewport(layout = grid.layout(2, 2), clip="off"))

      text.size <- 10

      vp <- vplayout(1,1)
      vp11 <- vp
      plot.2d <- plot.2d + theme(text = element_text(size = text.size))
      print(plot.2d, vp=vp)

      if(!(is.null(positionsToHighlight))) {
        # Merge the data and the positions to highlight by chr (col 1)
        # and start (col 2)
        addpts = merge(v, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))        
        frequencies.highlight <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, row.names=NULL, stringsAsFactors=NULL)

        # Add the labels
        if(length(addpts[,1]) > 1){
          if(highlightsHaveNames){

            # This code adapted from stackoverflow.com/questions/10536396/using-grconvertx-grconverty-in-ggplot2
            # Create a new viewport with clipping disabled so we can
            # put text outside the plot
            depth <- downViewport('panel.3-4-3-4')
            pushViewport(dataViewport(xData=c(0,100), yData=c(0,100), clip='off'))
            for(i in 1:dim(addpts)[1]) {
              if(addpts$vaf.1[i] < 1) {
                x <- addpts$vaf.1[i] - 12
                y <- addpts$vaf.2[i]
              } else if(addpts$vaf.2[i] < 1) {
                x <- addpts$vaf.1[i]
                y <- addpts$vaf.2[i] - 10                
              } else {
                x <- addpts$vaf.1[i] + 5
                y <- addpts$vaf.2[i] + 5
              }
              label <- addpts$gene_name[i]
              df <- data.frame(x=x, y=y)
              grid.text(x=x,y=y,label=label,default.units="native", gp=gpar(fontsize=8))
            }

            # Move up depth+1 in the tree.  NB: +1 because we pushed a
            # viewport on to the tree.
            upViewport(depth+1)
          }
        }        
        
      }

      vp <- vplayout(1,2)
      plot.1d.2 <- plot.1d.2 + theme(text = element_text(size = text.size))
      plot.1d.2 <- plot.1d.2 + coord_flip(ylim = c(xmin, res.1d.max.densities[[d2]]*1.1), xlim = c(xmin, xmax))
      print(plot.1d.2, vp=vp)

      vp <- vplayout(2,1)
      plot.1d.1 <- plot.1d.1 + theme(text = element_text(size = text.size))
      print(plot.1d.1, vp=vp)

      
    }
  }
  devoff = dev.off()
}



##---------------------------------------------------------------------
## create the two dimensional plot with scatter annotated with
## clustering results and 1D plots along margins
##---------------------------------------------------------------------------------
## Create two dimensional plot with scatter annotated with clustering result
##
sc.plot2d <- function(sco, outputFile, positionsToHighlight=NULL, highlightsHaveNames=FALSE, overlayClusters=TRUE, ellipse.metadata = list()){
  pdf(outputFile, width=7.2, height=6, bg="white")

  densityData = sco@densities
  vafs.merged = sco@vafs.merged
  sampleNames = sco@sampleNames
  dimensions = sco@dimensions
  clust = sco@clust

  numClusters = 0
  if(!is.null(vafs.merged$cluster)) {
    numClusters = max(vafs.merged$cluster, na.rm=T)
  }
  library(plotrix)


  ##create a 2d plot for each pairwise combination of samples
  for(d1 in 1:(dimensions-1)){
    for(d2 in d1:dimensions){
      if(d1==d2){
        next
      }

      vafs1 = getOneSampleVafs(vafs.merged, d1, numClusters);
      vafs2 = getOneSampleVafs(vafs.merged, d2, numClusters);

      ##get only cn2 points
      vafs1 = vafs1[vafs1$cleancn==2,]
      vafs2 = vafs2[vafs2$cleancn==2,]

      if(!is.null(vafs.merged$cluster)) {
        v = merge(vafs1,vafs2,by.x=c(1,2,8), by.y=c(1,2,8),suffixes=c(".1",".2"))
      } else {
        v = merge(vafs1,vafs2,by.x=c(1,2), by.y=c(1,2),suffixes=c(".1",".2"))
      }

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

      # If we will be highlighting some points, exclude them from
      # the general list of points to plot and plot them instead with
      # a different symbol/color (a black *)
      v.no.highlight <- v
      if(!(is.null(positionsToHighlight))) {
        chr.start.v <- cbind(v[,1], v[,2])
        chr.start.highlight <- cbind(positionsToHighlight[,1], positionsToHighlight[,2])        
        v.no.highlight <- v[!(apply(chr.start.v, 1, paste, collapse="$$") %in% apply(chr.start.highlight, 1, paste, collapse="$$")),]

      }
      
      if(!is.null(vafs.merged$cluster)) {      
        for(i in 1:numClusters){
          if(overlayClusters){
            if(dim(v.no.highlight[v.no.highlight$cluster==i,])[1] > 0) {
              points(v.no.highlight[v.no.highlight$cluster==i,]$vaf.1, v.no.highlight[v.no.highlight$cluster==i,]$vaf.2, col=cols[i], pch=i)
            }
          } else {
            if(dim(v.no.highlight[v.no.highlight$cluster==i,])[1] > 0) {
              points(v.no.highlight[v.no.highlight$cluster==i,]$vaf.1, v.no.highlight[v.no.highlight$cluster==i,]$vaf.2, pch=14)
            }
          }
        }
      } else {
        if(dim(v.no.highlight)[1] > 0) {
          points(v.no.highlight$vaf.1, v.no.highlight$vaf.2, pch=14)
        }
      }

      # Now plot the highlighted points so they are overlaid
      if(!(is.null(positionsToHighlight))) {
        # Merge the data and the positions to highlight by chr (col 1)
        # and start (col 2)
        addpts = merge(v, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))        
        
        # Plot the highlighted items.  NB:  we never overlay the
        # cluster on them, but expect this will be obvious from context
        points(addpts$vaf.1, addpts$vaf.2, pch="*", col="black", cex=2)

      }

      
      if(!is.null(vafs.merged$cluster)) {            
        for(i in 1:numClusters){
          if((!is.null(ellipse.metadata$SEMs.lb)) & (!is.null(ellipse.metadata$SEMs.ub))) {
            xc <- ellipse.metadata$SEMs.lb[i,d1] + ((ellipse.metadata$SEMs.ub[i,d1] - ellipse.metadata$SEMs.lb[i,d1])/2)
            yc <- ellipse.metadata$SEMs.lb[i,d2] + ((ellipse.metadata$SEMs.ub[i,d2] - ellipse.metadata$SEMs.lb[i,d2])/2)

            # ell <- my.ellipse(hlaxa = ((ellipse.metadata$SEMs.ub[i,d1] - ellipse.metadata$SEMs.lb[i,d1])/2), hlaxb = ((ellipse.metadata$SEMs.ub[i,d2] - ellipse.metadata$SEMs.lb[i,d2])/2), xc = xc, yc = yc)

            draw.ellipse(xc, yc, a = ((ellipse.metadata$SEMs.ub[i,d1] - ellipse.metadata$SEMs.lb[i,d1])/2), b = ((ellipse.metadata$SEMs.ub[i,d2] - ellipse.metadata$SEMs.lb[i,d2])/2))

          }

          if((!is.null(ellipse.metadata$std.dev.lb)) & (!is.null(ellipse.metadata$std.dev.ub))) {
            xc <- ellipse.metadata$std.dev.lb[i,d1] + ((ellipse.metadata$std.dev.ub[i,d1] - ellipse.metadata$std.dev.lb[i,d1])/2)
            yc <- ellipse.metadata$std.dev.lb[i,d2] + ((ellipse.metadata$std.dev.ub[i,d2] - ellipse.metadata$std.dev.lb[i,d2])/2)

            # Plot std dev as dashed line.
            draw.ellipse(xc, yc, a = ((ellipse.metadata$std.dev.ub[i,d1] - ellipse.metadata$std.dev.lb[i,d1])/2), b = ((ellipse.metadata$std.dev.ub[i,d2] - ellipse.metadata$std.dev.lb[i,d2])/2), lty=2)
          }
        }
      }
      #plot(-100, -100, xlim=c(0,100), ylim=c(0,100),main="Clusters")
      if(!is.null(vafs.merged$cluster)) {
        legend("topright", legend=1:numClusters, col=cols[1:numClusters], title="Clusters", pch=1:numClusters)
      }

      # Add annotation for gene names, if requested
      if(!(is.null(positionsToHighlight))){
        addpts = merge(v, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))
        # write.table(file="genes.txt", unique(addpts$gene_name), row.names=FALSE, col.names=FALSE, quote=FALSE)
        if(length(addpts[,1]) > 1){
          if(highlightsHaveNames){
            for(i in 1:dim(addpts)[1]) {
              par(xpd=NA)
              cex <- 1
              if(addpts$vaf.1[i] < 1) {
                text(addpts$vaf.1[i] - 8,addpts$vaf.2[i],labels=addpts$gene_name[i],cex=cex)
              } else if(addpts$vaf.2[i] < 1) {
                text(addpts$vaf.1[i],addpts$vaf.2[i] - 5,labels=addpts$gene_name[i],cex=cex)
              } else { 
                text(addpts$vaf.1[i] + 4,addpts$vaf.2[i] + 4,labels=addpts$gene_name[i],cex=cex)
              }
            }
          }
        }
      } # End add gene annotations
    }
  }
  devoff = dev.off()
}


##-------------------------------------------------------------------------------------
## plot three samples in 3d, optionally create a GIF
##
sc.plot3d <- function(sco, samplesToPlot, size=700, outputFile=NULL){
  library(rgl)
  ##set the size of the window
  r3dDefaults$windowRect <- c(0,50, size, size)

  if(length(samplesToPlot) != 3){
    print("ERROR: must provide exactly 3 sample names for 3d plotting")
    return(0)
  }

  a = sco@vafs.merged[,c(paste(samplesToPlot,".vaf",sep=""),"cluster")]
  a = a[!is.na(a$cluster),]
  numClusters=length(unique(a$cluster))
  cols=getClusterColors(numClusters)
  colvec = cols[a$cluster]

  plot3d(a[,1], a[,2], a[,3], xlim=c(0,100), ylim=c(0,100),zlim=c(0,100), axes=FALSE,
         xlab=samplesToPlot[1], ylab=samplesToPlot[2], zlab=samplesToPlot[3],
         type="s", col=colvec)
  ##add a box
  axes3d( edges=c("x--", "y--", "z"),labels=FALSE)
  for(i in c("+","-")){
    for(j in c("+","-")){
      axes3d( edges=paste("x",i,j,sep=""), tick=FALSE, labels=FALSE)
      axes3d( edges=paste("y",i,j,sep=""), tick=FALSE, labels=FALSE)
      axes3d( edges=paste("z",i,j,sep=""), tick=FALSE, labels=FALSE)
    }
  }

  if(is.null(outputFile)){
    play3d(spin3d(axis=c(0,0,1), rpm=10), duration=6)
  } else {
    ##remove trailing .gif, since movie3d adds it
    outputFile = sub(".gif$","",outputFile)
    movie3d(spin3d(axis=c(0,0,1), rpm=10), duration=6, dir=getwd(), movie=outputFile)
  }
}


##-------------------------------------------------------------------------------------
## Get a list of colors to use for the clusters (after 20 they start repeating)
##
getClusterColors <- function(numClusters){
  library(RColorBrewer)
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
  if(numClusters > 0){
    vafs = cbind(vafs, vafs.merged$cluster)
    names(vafs)[8] = "cluster"
  }
  return(vafs)
}



