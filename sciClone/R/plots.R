#---------------------------------------------------------------------------------
## Create the one dimensional plot with kde and scatter
##
sc.plot1d <- function(sco, outputFile,
                   plotOnlyCN2=FALSE, showCopyNumberScatterPlots=TRUE, highlightSexChrs=TRUE,
                   positionsToHighlight=NULL, highlightsHaveNames=FALSE, overlayClusters=TRUE,
                   overlayIndividualModels=TRUE, showHistogram=FALSE, onlyLabelHighestPeak=FALSE,
                   minimumLabelledPeakHeight=0.001, showTitle=TRUE){


  densityData = sco@densities
  vafs.merged = sco@vafs.merged

  sampleNames = sco@sampleNames
  dimensions = sco@dimensions

  ##are we plotting everything or just CN2?
  cnToPlot = c();
  if(plotOnlyCN2){
    cnToPlot = c(2)
  } else {
    cnToPlot = 1:4
  }

  # If any of the vafs are named, assume we will be plotting them and
  # will need a legend for them.
  if(!is.null(positionsToHighlight)) {
    names(positionsToHighlight)=c("chr","st","name");
  }
  addpts = merge(vafs.merged, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))
  add.legend <- FALSE
  if((dim(addpts)[1] > 0) & (any(addpts$name != ""))) { 
    if(showCopyNumberScatterPlots & ( length(cnToPlot) < 2 ) & highlightsHaveNames) {
      add.legend <- TRUE
      # We will add a legend panel
    }
  }

  #sanity checks
  if(highlightsHaveNames & add.legend){
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

  clust = NULL
  if(overlayClusters){
    if(is.null(sco@clust[1])){
      print("ERROR: can't overlay clusters when clustering was not done on the input data")
      return(0)
    } else {
      clust = sco@clust
    }
  }

  num.rows <- length(cnToPlot) + 1

  if(add.legend) {
    num.rows <- num.rows + 1
  }
  
  # 3.3 x 7.5 is a good dimensionality for 5 rows.  Scale accordingly
  # if we have fewer rows.
  #pdf(file=outputFile, width=3.3, height=7.5, bg="white");
  height <- 7.5
  scale <- 5 / num.rows
  width <- 3.3 * scale
  pdf(file=outputFile, width=width, height=height, bg="white");

  numClusters = 0
  if(!is.null(clust)) {
    numClusters = max(clust$cluster.assignments)
  }

  ##one plot for each sample
  for(d in 1:dimensions){

    name=sampleNames[d]

    ##grab only the vafs for this sample:
    vafs = getOneSampleVafs(vafs.merged, name, numClusters);

    par(mfcol=c(num.rows,1),mar=c(0.5*scale,3*scale,1*scale,1.5*scale),oma=c(3,0,4,0),mgp = c(3,1,0));
    
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
    axis(side=2,at=c(0,1),labels=c(0, ""), las=1, cex.axis=0.6*scale, hadj=0.6,
         lwd=0.5*scale, lwd.ticks=0.5*scale, tck=-0.01);


    ##colors for different copy numbers
    colors=c("#1C3660AA","#67B32EAA","#F49819AA","#E52420AA")

    density.curve.width <- 4
    for(i in cnToPlot){
      if(!(is.null(densities[[i]])) & (!showHistogram | (i!= 2))){
        ##density lines
        lines(densities[[i]]$x, scalingFactor*factors[[i]], col=colors[i], lwd=density.curve.width*scale);
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
                 cex=0.7*scale, srt=0, col=colors[[i]]);
          }
        }
      } else if(showHistogram & (i == 2)) {

        ## Only show histogram for copy number = 2
        v = vafs[which(vafs$cleancn==2 & vafs$adequateDepth==1),];

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
    model.width <- density.curve.width / 2
    ## Plot the individual models with dotted lines (3) or dashed (2)
    individual.model.style <- 3
    individual.model.width <- density.curve.width / 2

    if(!(is.null(clust))){
      maxFitDensity <- max(clust$fit.y[d,])
      #points(clust$fit.x, clust$fit.y[d,]*25, type="l",col="grey50")
      lines(clust$fit.x, clust$fit.y[d,]/maxFitDensity, type="l",col="grey50",lty=model.style, lwd=model.width*scale)
      if(overlayIndividualModels==TRUE) {
        for(i in 1:numClusters) {
          lines(clust$fit.x, clust$individual.fits.y[[i]][d,]/maxFitDensity,
                type="l",col="grey50",lty=individual.model.style, lwd=individual.model.width*scale)
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
      if(showHistogram == FALSE) {
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
    legend(x="topright", lwd=lwd, lty=lty, legend=leg, col=lcol, bty="n", cex=0.6*scale, y.intersp=1.25, pch=pchs, pt.bg = pt.bgs);


    axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6*scale,lwd=0.5*scale,lwd.ticks=0.5*scale,padj=1.4);
    mtext("Variant Allele Frequency",adj=0.5,padj=-3.1/scale,cex=0.6*scale,side=3);
    mtext("Density (a.u.)",side=2,cex=0.6*scale,padj=-4.2/scale);


    ##add a title to the plot
    if(showTitle){
      title=""
      if(is.null(sampleNames[d])){
        title="Clonality Plot"
      } else {
        title=paste(sampleNames[d],"Clonality Plot",sep=" ");
      }
      mtext(title, adj=0.5, padj=-5/scale, cex=0.65*scale, side=3);
    }

    ##-----------------------------------------------------
    ##create the scatterplots of vaf vs copy number

    if(showCopyNumberScatterPlots) {
      for(i in cnToPlot){
        v = vafs[which(vafs$cleancn==i & vafs$adequateDepth==1),];
        drawScatterPlot(v, highlightSexChrs, positionsToHighlight, colors, i, maxDepth, highlightsHaveNames, overlayClusters, scale)
        axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6*scale,lwd=0.5*scale,lwd.ticks=0.5*scale,padj=-1.4);


        if(length(cnToPlot) < 2 & highlightsHaveNames){
          addHighlightLegend(v, positionsToHighlight,scale)
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
drawScatterPlot <- function(data, highlightSexChrs, positionsToHighlight, colors, cn, maxDepth, highlightsHaveNames, overlayClusters,scale=1){

  ## define plot colors
  ptcolor = colors[cn]
  circlecolor = substr(colors[cn],1,7) #chop off the alpha value

  ## define the plot space by plotting offscreen points
  plot.default( x=-10000, y=1, log="y", type="p", pch=19, cex=0.4*scale,
               col="#00000000", xlim=c(0,100), ylim=c(5,maxDepth*3),
               axes=FALSE, ann=FALSE, xaxs="i", yaxs="i");

  addPoints <- function(data, color, highlightSexChrs, pch=NULL){
    outlineCol = rgb(0,0,0,0.1);
    if(highlightSexChrs){
      ##plot autosomes
      data.autosomes = data[!(data$chr == "X" | data$chr == "Y"),]
      points(data.autosomes$vaf, data.autosomes$depth, type="p", pch=16, cex=0.75*scale, col=color);
      #points(data.autosomes$vaf, data.autosomes$depth, type="p", pch=1, cex=0.75*scale, col=outlineCol, lwd=);
      ##plot sex chromsomes with different shape
      data.sex = data[(data$chr == "X" | data$chr == "Y"),]
      points(data.sex$vaf, data.sex$depth, type="p", pch=17, cex=0.75*scale, col=color);
      #points(data.sex$vaf, data.sex$depth, type="p", pch=1, cex=0.75*scale, col=outlineCol);
    } else {
      points(data$vaf, data$depth, type="p", pch=16, cex=0.75*scale, col=color);
      ##add outline
      #points(data$vaf, data$depth, type="p", pch=1, cex=0.75*scale, col=outlineCol);
    }

  }

  #do we have any points to plot?
  if(length(data[,1]) > 0){
    ##if we have cluster assignments in col 8, color them
    if(any(grepl(pattern="^cluster$",names(data))) & overlayClusters & cn==2){
      numClusters=max(data[,c("cluster")],na.rm=T)
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
      addpts = merge(data, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))
      if(length(addpts[,1]) > 1){
        if(highlightsHaveNames){
          for(i in 1:length(addpts$vaf)){
            text(addpts$vaf[i],addpts$depth[i],labels=i,cex=0.5*scale)
          }
        } else {
          addPoints(addpts, col="#555555FF", highlightSexChrs);
          ##points(x=addpts$vaf,y=addpts$depth,type="p",pch=7,cex=0.8*scale,col="#555555FF");
        }
      }
    }
  }

  ## define the axis
  axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6*scale,hadj=0.5);
  for (i in 2:length(axTicks(2)-1)) {
    lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
  }

  ## plot the background color
  rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA);

  ## add cn circle
  points(x=c(97),y=c(maxDepth),type="p",pch=19,cex=3*scale,col=circlecolor);
  text(c(97),y=c(maxDepth), labels=c(cn), cex=1*scale, col="#FFFFFFFF")


  ## y axis label
  mtext("Tumor Coverage",side=2,cex=0.6*scale,padj=-4.2/scale);
}


##--------------------------------------------------------------------
## add a legend for highlighted points with names
addHighlightLegend <- function(data, positionsToHighlight, scale){

  if((is.null(positionsToHighlight))){ return() }

  names(positionsToHighlight)=c("chr","st","name");

  addpts = merge(data, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))

  if(length(addpts[,1]) == 0){ return() }

  non.trivial.names <- addpts$name[addpts$name != ""]

  if(length(non.trivial.names) == 0){ return() }
  
  plot.default(x=-10000, y=1, type="p", pch=19, cex=0.4,
               col="#00000000", xlim=c(0,1000), ylim=c(0,1000),
               axes=FALSE, ann=FALSE);

  # mtext("Genes",side=2,cex=0.5,padj=-4.2);

  ypos=rev(seq(0,900,(900/13)))[1:13]
  # Don't put empty names on legend
  non.trivial.indices <- (1:length(addpts$name))[addpts$name != ""]
  ncol=ceiling(length(non.trivial.names)/13)
  xpos=0;
  offset=1
  nxt <- 1
  for(n in 1:ncol){
    names = non.trivial.names[offset:(offset+12)]
    names = as.character(names[!(is.na(names))])
    num = length(names)

    for(i in 1:num){
      text(xpos, ypos[i], paste(non.trivial.indices[nxt],". ",names[i],sep=""), cex=0.6*scale, pos=4)
      nxt <- nxt+1
    }
    xpos=xpos+250;
    offset=offset+13;
  }
}


##---------------------------------------------------------------------
## create the two dimensional plot with scatter annotated with
## clustering results and 1D plots along margins, this time using
## ggplot2

plot2dWithMargins <- function(sco, outputFile,positionsToHighlight=NULL, highlightsHaveNames=FALSE, overlayErrorBars=FALSE) {
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

    vafs = getOneSampleVafs(vafs.1d.merged[[d]], sampleNames[d], numClusters)

    # Only show copy number = 2
    v = vafs[which(vafs$cleancn==2 & vafs$adequateDepth==1),];

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

    # If we actually performed clustering, show it.
    if(numClusters > 0) {
      f <- splinefun(marginalClust[[d]]$fit.x, scale*marginalClust[[d]]$fit.y[1,,drop=FALSE])

      g <- g + stat_function(data = limits, fun=f, mapping=aes(x))
    }
      
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

      vafs1 = getOneSampleVafs(vafs.merged, sampleNames[d1], numClusters);
      vafs2 = getOneSampleVafs(vafs.merged, sampleNames[d2], numClusters);

      ##get only cn2 points
      vafs1 = vafs1[which(vafs1$cleancn==2 & vafs1$adequateDepth==1),]
      vafs2 = vafs2[which(vafs2$cleancn==2 & vafs2$adequateDepth==1),]

      if(numClusters > 0) {
        v = merge(vafs1,vafs2,by.x=c("chr","st","cluster"), by.y=c("chr","st","cluster"),suffixes=c(".1",".2"))
        # Remove any outliers--these will have cluster assignment 0
        v.outlier <- v[v$cluster == 0,]
        v <- v[v$cluster != 0,]
      } else {
        v = merge(vafs1,vafs2,by.x=c("chr","st"), by.y=c("chr","st"),suffixes=c(".1",".2"))
      }
      
      v.no.highlight <- v
      if(!(is.null(positionsToHighlight))) {
        names(positionsToHighlight)=c("chr","st","name");
        chr.start.v <- cbind(v[,"chr"], v[,"st"])
        chr.start.highlight <- cbind(positionsToHighlight[,1], positionsToHighlight[,2])
        v.no.highlight <- v[!(apply(chr.start.v, 1, paste, collapse="$$") %in% apply(chr.start.highlight, 1, paste, collapse="$$")),]

      }

      title <- ""


      # Plot the points that we will not highlight
      frequencies.no.highlight <- data.frame(x=v.no.highlight$vaf.1, y=v.no.highlight$vaf.2, row.names=NULL, stringsAsFactors=NULL)

      if(numClusters > 0) {
        clusters <- v.no.highlight$cluster
        cols=getClusterColors(numClusters)
        colvec = cols[clusters]

        g <- ggplot(data = frequencies.no.highlight, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = frequencies.no.highlight, aes(x=x, y=y), shape=clusters, colour=colvec)
      } else {
        g <- ggplot(data = frequencies.no.highlight, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = frequencies.no.highlight, aes(x=x, y=y))
      }
      
      if(overlayErrorBars == TRUE) {
        err.bars.1 <- compute.binomial.error.bars(v.no.highlight$var.1, v.no.highlight$depth.1) * 100
        err.bars.2 <- compute.binomial.error.bars(v.no.highlight$var.2, v.no.highlight$depth.2) * 100
        err.df.x <- data.frame(x=v.no.highlight$vaf.1, y=v.no.highlight$vaf.2, xmin=err.bars.1$lb, xmax=err.bars.1$ub)
        g <- g + geom_errorbarh(data = err.df.x, aes(x=x, y=y, xmin=xmin, xmax=xmax), colour=colvec)

        err.df.y <- data.frame(x=v.no.highlight$vaf.1, y=v.no.highlight$vaf.2, ymin=err.bars.2$lb, ymax=err.bars.2$ub)
        g <- g + geom_errorbar(data = err.df.y, aes(x=x, y=y, ymin=ymin, ymax=ymax), colour=colvec)
      } 
      

      # Now overlay any points that we will highlight
      if(!(is.null(positionsToHighlight))) {
        # Merge the data and the positions to highlight by chr (col 1)
        # and start (col 2)
        addpts = merge(v, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))
        if(dim(addpts)[1] > 0) {
          frequencies.highlight <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, row.names=NULL, stringsAsFactors=NULL)

          g <- g + geom_point(data = frequencies.highlight, aes(x=x, y=y), shape="*", size=10, colour="black")

          if(overlayErrorBars == TRUE) {
            err.bars.1 <- compute.binomial.error.bars(addpts$var.1, addpts$depth.1) * 100
            err.bars.2 <- compute.binomial.error.bars(addpts$var.2, addpts$depth.2) * 100
            err.df.x <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, xmin=err.bars.1$lb, xmax=err.bars.1$ub)
            g <- g + geom_errorbarh(data = err.df.x, aes(x=x, y=y, xmin=xmin, xmax=xmax), colour="black")

            err.df.y <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, ymin=err.bars.2$lb, ymax=err.bars.2$ub)
            g <- g + geom_errorbar(data = err.df.y, aes(x=x, y=y, ymin=ymin, ymax=ymax), colour="black")
          }
        } 

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
        addpts = merge(v, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))
        frequencies.highlight <- data.frame(x=addpts$vaf.1, y=addpts$vaf.2, row.names=NULL, stringsAsFactors=NULL)

        # Add the labels
        if(dim(addpts)[1] > 0){
          if(highlightsHaveNames){

            # This code adapted from stackoverflow.com/questions/10536396/using-grconvertx-grconverty-in-ggplot2
            # Create a new viewport with clipping disabled so we can
            # put text outside the plot
            depth <- downViewport('panel.3-4-3-4')
            pushViewport(dataViewport(xData=c(0,100), yData=c(0,100), clip='off'))

            xs <- list()
            ys <- list()
            labels <- list()

            nxt <- 1
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
              label <- as.character(addpts$name[i])
              #df <- data.frame(x=x, y=y)
              if(label == "") {
                # No label to place
              } else { 
                # grid.text(x=x,y=y,label=label,default.units="native", gp=gpar(fontsize=8))
                xs[[nxt]] <- x
                ys[[nxt]] <- y
                labels[[nxt]] <- label
                nxt <- nxt + 1
              }
            }

            xs <- unlist(xs)
            ys <- unlist(ys)
            labels <- unlist(labels)
            num.labels <- length(labels)
            library(TeachingDemos)
            # By using the min argument, ensure that the annotations
            # do not overlap the corresponding symbol.  NB:  we anticipate
            # this code will only be active for the interior points
            if(num.labels > 0){             
              xs <- spread.labs(xs, mindiff=4, min=xs)
              ys <- spread.labs(ys, mindiff=4, min=ys)
              for(i in 1:num.labels) {
                grid.text(x=xs[i],y=ys[i],label=labels[i],default.units="native", gp=gpar(fontsize=8))
              }
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


## Compute the lower and upper bound for a "1-std dev" binomial confidence
## interval for a given number of successes and total number trials
compute.binomial.error.bars <- function(successes, total.trials){
  suppressPackageStartupMessages(library(MKmisc))
  suppressPackageStartupMessages(library(NORMT3))
  # Return a "1 std dev" confidence interval
  width <- as.real(erf(1/sqrt(2)))
  lb <- mapply(function(a,b) binomCI(a, b, conf.level=width, method="jeffreys")$CI[1], successes, total.trials)
  ub <- mapply(function(a,b) binomCI(a, b, conf.level=width, method="jeffreys")$CI[2], successes, total.trials)
  res <- data.frame(lb=lb, ub=ub)
  return(res)
}

##---------------------------------------------------------------------
## create the two dimensional plot with scatter annotated with
## clustering results and 1D plots along margins
##---------------------------------------------------------------------------------
## Create two dimensional plot with scatter annotated with clustering result
##
sc.plot2d <- function(sco, outputFile, positionsToHighlight=NULL, highlightsHaveNames=FALSE, overlayClusters=TRUE, overlayErrorBars=FALSE, ellipse.metadata = list()){
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
  library(plotrix)  # For plotCI among others.


  ##create a 2d plot for each pairwise combination of samples
  for(d1 in 1:(dimensions-1)){
    for(d2 in d1:dimensions){
      if(d1==d2){
        next
      }
      vafs1 = getOneSampleVafs(vafs.merged, sampleNames[d1], numClusters);
      vafs2 = getOneSampleVafs(vafs.merged, sampleNames[d2], numClusters);
      ##get only cn2 points with adequate coverage
      vafs1 = vafs1[which(vafs1$cleancn==2 & vafs1$adequateDepth==1),]
      vafs2 = vafs2[which(vafs2$cleancn==2 & vafs2$adequateDepth==1),]

      if(!is.null(vafs.merged$cluster)) {
        v = merge(vafs1,vafs2,by.x=c("chr","st","cluster"), by.y=c("chr","st","cluster"),suffixes=c(".1",".2"))
        # Remove any outliers--these will have cluster assignment 0
        v.outlier <- v[v$cluster == 0,]
        v <- v[v$cluster != 0,]
      } else {
        v = merge(vafs1,vafs2,by.x=c("chr","st"), by.y=c("chr","st"),suffixes=c(".1",".2"))
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
        names(positionsToHighlight)=c("chr","st","name");        
        chr.start.v <- cbind(v[,"chr"], v[,"st"])
        chr.start.highlight <- cbind(positionsToHighlight[,1], positionsToHighlight[,2])
        v.no.highlight <- v[!(apply(chr.start.v, 1, paste, collapse="$$") %in% apply(chr.start.highlight, 1, paste, collapse="$$")),]

      }

      if(overlayErrorBars == TRUE) {
        err.bars.1 <- compute.binomial.error.bars(v.no.highlight$var.1, v.no.highlight$depth.1) * 100
        err.bars.2 <- compute.binomial.error.bars(v.no.highlight$var.2, v.no.highlight$depth.2) * 100
      } 

      if(!is.null(vafs.merged$cluster)) {      
        for(i in 1:numClusters){
          indices <- v.no.highlight$cluster==i
          if(overlayClusters){
            if(dim(v.no.highlight[indices,])[1] > 0) {
              if(overlayErrorBars == TRUE) { 
                plotCI(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, col=cols[i], pch=i, li=err.bars.1$lb[indices], ui=err.bars.1$ub[indices], add=TRUE, err="x")
                plotCI(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, col=cols[i], pch=i, li=err.bars.2$lb[indices], ui=err.bars.2$ub[indices], add=TRUE, err="y")
              } else {
                points(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, col=cols[i], pch=i)
              }
            }
          } else {
            if(dim(v.no.highlight[indices,])[1] > 0) {
              if(overlayErrorBars == TRUE) { 
                plotCI(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, pch=14, li=err.bars.1$lb[indices], ui=err.bars.1$ub[indices], add=TRUE, err="x")
                plotCI(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, pch=14, li=err.bars.2$lb[indices], ui=err.bars.2$ub[indices], add=TRUE, err="y")
              } else {
                points(v.no.highlight[indices,]$vaf.1, v.no.highlight[indices,]$vaf.2, pch=14)
              }
            }
          }
        }
      } else {
        if(dim(v.no.highlight)[1] > 0) {
          if(overlayErrorBars == TRUE) { 
            plotCI(v.no.highlight$vaf.1, v.no.highlight$vaf.2, pch=14, li=err.bars.1$lb, ui=err.bars.1$ub, add=TRUE, err="x")
            plotCI(v.no.highlight$vaf.1, v.no.highlight$vaf.2, pch=14, li=err.bars.2$lb, ui=err.bars.2$ub, add=TRUE, err="y")
          } else {
            points(v.no.highlight$vaf.1, v.no.highlight$vaf.2, pch=14)
          }
        }
      }

      # Now plot the highlighted points so they are overlaid
      if(!(is.null(positionsToHighlight))) {
        # Merge the data and the positions to highlight by chr (col 1)
        # and start (col 2)
        addpts = merge(v, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))

        # Plot the highlighted items.  NB:  we never overlay the
        # cluster on them, but expect this will be obvious from context
        if(dim(addpts)[1] > 0) {
          if(overlayErrorBars == TRUE) { 
            err.bars.1 <- compute.binomial.error.bars(addpts$var.1, addpts$depth.1) * 100
            err.bars.2 <- compute.binomial.error.bars(addpts$var.2, addpts$depth.2) * 100
            plotCI(addpts$vaf.1, addpts$vaf.2, pch="*", col="black", cex=2, li=err.bars.1$lb, ui=err.bars.1$ub, add=TRUE, err="x")
            plotCI(addpts$vaf.1, addpts$vaf.2, pch="*", col="black", cex=2, li=err.bars.2$lb, ui=err.bars.2$ub, add=TRUE, err="y")
          } else {
            points(addpts$vaf.1, addpts$vaf.2, pch="*", col="black", cex=2)
          }
        }

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
        addpts = merge(v, positionsToHighlight, by.x=c("chr","st"), by.y = c("chr","st"))
        # write.table(file="genes.txt", unique(addpts$gene_name), row.names=FALSE, col.names=FALSE, quote=FALSE)
        if(dim(addpts)[1] > 0){
          if(highlightsHaveNames){
            xs <- list()
            ys <- list()
            labels <- list()

            nxt <- 1
            for(i in 1:dim(addpts)[1]) {
              par(xpd=NA)
              cex <- 1
            
              if(addpts$vaf.1[i] < 1) {
                #text(addpts$vaf.1[i] - 8,addpts$vaf.2[i],labels=addpts$gene_name[i],cex=cex)
                x <- addpts$vaf.1[i] - 8
                y <- addpts$vaf.2[i]
              } else if(addpts$vaf.2[i] < 1) {
                #text(addpts$vaf.1[i],addpts$vaf.2[i] - 5,labels=addpts$gene_name[i],cex=cex)
                x <- addpts$vaf.1[i]
                y <- addpts$vaf.2[i] - 5                
              } else { 
                #text(addpts$vaf.1[i] + 4,addpts$vaf.2[i] + 4,labels=addpts$gene_name[i],cex=cex)
                x <- addpts$vaf.1[i] + 4
                y <- addpts$vaf.2[i] + 4
              }
              label <- as.character(addpts$name[i])
              #df <- data.frame(x=x, y=y)
              #grid.text(x=x,y=y,label=label,default.units="native", gp=gpar(fontsize=8))
              if(label != "") {
                xs[[nxt]] <- x
                ys[[nxt]] <- y
                labels[[nxt]] <- label
                nxt <- nxt + 1
              }
            }

            xs <- unlist(xs)
            ys <- unlist(ys)
            labels <- unlist(labels)
            num.labels <- length(labels)
            library(TeachingDemos)
            # By using the min argument, ensure that the annotations
            # do not overlap the corresponding symbol.  NB:  we anticipate
            # this code will only be active for the interior points
            if(num.labels > 0){ 
              xs <- spread.labs(xs, mindiff=4, min=xs)
              ys <- spread.labs(ys, mindiff=4, min=ys)
              for(i in 1:num.labels) {
                #grid.text(x=xs[i],y=ys[i],label=labels[i],default.units="native", gp=gpar(fontsize=8))
                text(x=xs[i], y=ys[i], label=labels[i], cex=cex)
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
  numClusters=max(a$cluster)
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
getOneSampleVafs <- function(vafs.merged, name, numClusters){
  common = c("chr","st","adequateDepth")
  a = vafs.merged[,common]
  prefix = paste("^",name,".",sep="")
  cols=grepl(prefix, names(vafs.merged))
  header=sub(prefix, "", names(vafs.merged)[cols])
  vafs = cbind(a,vafs.merged[,cols])
  names(vafs) = c(common,header)

  if(numClusters > 0){
    vafs = cbind(vafs, vafs.merged$cluster)
    names(vafs)[length(vafs)] = "cluster"
  }
  return(vafs)
}



