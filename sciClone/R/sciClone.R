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
sciClone <- function(vafs, outputImage, copyNumberCalls=NULL, sampleName,
                        minimumDepth=100, clusteredDataOutputFile=NULL,
                        componentDistribution="Binomial", maximumClustersToTest=8,
                        minimumLabelledPeakHeight=0.001, onlyLabelHighestPeak=FALSE,
                        overlayClusters=FALSE, plotOnlyCN2=FALSE, positionsToHighlight=NULL,
                        purity=0, highlightSexChrs=TRUE, testing=FALSE, cnCallsAreLog2=FALSE,
                        useSexChrs=TRUE, highlightsHaveNames=FALSE, doClustering=TRUE){

  library(mixtools)
  library(mixdist)
  names(vafs) = c("chr","st","ref","var","vaf")

  ##remove MT values
  vafs = vafs[!(vafs$chr == "M" | vafs$chr == "MT"),]

  ##remove NA sites
  vafs = vafs[!(is.na(vafs$vaf)),]
  
  if(is.null(copyNumberCalls)){
    ##assume all sites are 2x if no cn info
    print("No copy number file specified. Assuming all variants have a CN of 2.")
    vafs$cn = 2;
  } else {
    if(cnCallsAreLog2){
      copyNumberCalls[,4] = (2^(copyNumberCalls[,4]))*2
    }
    vafs = addCnToVafs(vafs,copyNumberCalls)
  }

  ##sanity check
  if(highlightsHaveNames){
    if(is.null(positionsToHighlight)){
      print("ERROR - if highlightsHaveNames is true, positionsToHighlight must be provided")
      stop()
    }
    plotOnlyCN2=TRUE;
  }
  
  ##remove sex chromosomes if specified
  if(!(useSexChrs)){
    vafs = vafs[vafs$chr != "X" & vafs$chr != "Y",]
  }
  
  ##add depth
  vafs = vafs[vafs$vaf > 0,]
  vafs$depth = round(vafs$var/(vafs$vaf/100))

  ##remove any sites with less than the minimum depth
  vafs = vafs[vafs$depth >= minimumDepth,]
  if(length(vafs$chr) == 0){
    print(paste("No variants exceed a depth of",minimumDepth,". Lower this threshold and try again."))
    stop()
  }
  print(paste("Number of variants with depth >= ",minimumDepth," being used for analysis: ",length(vafs$chr),sep=""))

  
  
  ##calculate the densities and peaks for variants of each copy number
  ##default cutoffs are +/- 0.5x cn
  densities = vector("list",4)
  factors = vector("list",4)
  peakPos = vector("list",4)
  peakHeights = vector("list",4)
  vafsByCn = vector("list",4)

  cnLoThresh = c(0,1.5,2.5,3.5)
  cnHiThresh = c(1.5,2.5,3.5,4.5)

  maxDensity = 0
  maxDepth = 0;

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

      ##store the vaf info
      vafsByCn[[i]] = v

      ##store the largest depth for use in scaling the plots later
      if(max(vafsByCn[[i]]$depth) > maxDepth){
        maxDepth = max(vafsByCn[[i]]$depth)
      }

    } #else has a value of NULL
  }

  
  ## This stop should probably be replaced so that plotting can take place
  ## for tumors with ployplody, even if we can't cluster
  ## (maybe we should even cluster with 3x regions, etc - put it on the todo list)
  if(is.null(densities[[2]])){
    print("can't do clustering - no copy-number neutral regions to operate on");
    stop();
  }


  ##determine tumor purity if not specified
  if (purity == 0){    
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
    #if all of these failed to find a good peak, assume 100% purity
    if (purity == 0) {
      print("Unable to make reliable calculation of purity, so assuming 100%")
      purity = 100;
    }
    print(paste("Tumor purity estimated to be ",signif(purity,4),".",sep=""));
  }


  if(doClustering){
    ##----------------------------------------------------
    ## determine number of clusters in the dataset
    num_clusters = 0;
    print("Performing 'mixdist' analysis...");

    ##function to process mixdist results
    process_percents <- function(percents,chisq,pval,componentDistribution) {
      minchisq = NULL;
      best_fit_component_assessment = 0;
      minpval = NULL;
      true_cluster_count = NULL;

      for (i in 1:maximumClustersToTest) {
        if (percents[i]!="Error" && !is.nan(pval[i]) && pval[i] < 0.05) {
          if (is.null(minchisq) || chisq[i] < minchisq) {
            minchisq = chisq[i];
            minpval = pval[i];
            best_fit_component_assessment = i;
          }
        }
        true_cluster_count[i] = 0;
        percentage_per_cluster = as.numeric(percents[[i]]);
        true_clusters = percentage_per_cluster > 0.02;
        for (j in 1:i) {
          if (isTRUE(true_clusters[j])) {
            true_cluster_count[i] = true_cluster_count[i] + 1;
          }
        }
      }
      print(paste("chosen_component_assessment = ",best_fit_component_assessment,", true_component_count = ",true_cluster_count[best_fit_component_assessment],", assessment_chisq_value = ",minchisq,", assessment_p-value = ",minpval,sep=""));

      ## I was getting errors when no clusters exceeded the pval of 0.05
      ## in that case, it returns numeric(0), which causes problems later on
      ## solution is to return just 0 clusters and warn user
      if(best_fit_component_assessment > 0){
        return(true_cluster_count[best_fit_component_assessment]);
      } else {
        print("WARNING: unable to estimate number of clusters - significance of 0.05 not reached")
        return(0)
      }
    }
    

    ## run mixdist
    data = vafsByCn[[2]]$vaf
    grouped_data = NULL;

    ## if ceiling of max(data) is odd, then add 1 and group data. else, group data using the even ceiling
    if (ceiling(max(data))%%2) {
      grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data))+1,2)));
    } else {
      grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data)),2)));
    }

    ## for each component count, get mixdist fit estimates for normal distribution
    percents=NULL; chisq=NULL; pval=NULL;
    for (i in 1:maximumClustersToTest) {
      data_params = NULL;
      test = NULL;

      if (componentDistribution == "Normal") {
        data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data),i));
        test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="norm"), silent=TRUE);
      }
      if (componentDistribution == "Binomial") {
        data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data)/2,i));
        test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="binom", constr=mixconstr(consigma="BINOM",size=rep(round(length(data)),i))), silent=TRUE);
      }

      if (class(test) == 'try-error') {
        percents[i] = "Error";
        print(paste("Common mixdist error when looking for ",i," components.",sep=""));
      }
      else {
        percents[i]=list(test$parameters$pi)
        chisq[i]=test$chisq;
        pval[i] = test$P;

        if(testing){
          filename = paste("plot_component_",i,".pdf",sep="");
          print(paste("Testing output: plotting ",filename));
          pdf(filename)
          ##dev.new();
          plot(test);
          ##dev.copy(pdf,filename);
          d = dev.off();
        }
      }
    }
    num_clusters = process_percents(percents,chisq,pval,componentDistribution);
  } else {
    num_clusters = 0
  }

  

  ##--------------------------------------------------
  ##do the plotting
  pdf(file=outputImage,width=3.3,height=7.5,bg="white");
  par(mfcol=c(5,1),mar=c(0.5,3,1,1.5),oma=c(3,0,4,0),mgp = c(3,1,0));

  scalingFactor = 25/maxDensity;

  plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");
  ##plot bg color
  rect(0, 0, 100, 28, col = "#00000011",border=NA);

  axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxDensity)),las=1,cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);

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

  ##legend
  if( length(cnToPlot)== 1 ){
    legend(x="topright", lwd=2, legend=c("2 Copies"), col=colors[2], bty="n",
           cex=0.6, y.intersp=1.25);
  } else {
    legend(x="topright", lwd=2, legend=c("1 Copy","2 Copies","3 Copies","4 Copies"),
           col=colors, bty="n", cex=0.6, y.intersp=1.25);
  }


  ##legend for tier1 labeling
  
  

  axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
  mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
  mtext("Kernel Density",side=2,cex=0.5,padj=-4.2);

  ##add a title to the plot
  title=""
  if(is.null(sampleName)){
    title="Clonality Plot"
  } else {
    title=paste(sampleName,"Clonality Plot",sep=" ");
  }
  mtext(title, adj=0.5, padj=-5, cex=0.65, side=3);


  ##create the scatterplots of vaf vs density
  for(i in cnToPlot){

    drawScatterPlot(vafsByCn[[i]], highlightSexChrs, positionsToHighlight, colors, i, maxDepth, highlightsHaveNames)

    ## Plot clustered points
    if(i==2 && overlayClusters && num_clusters > 0){
      ## cluster tumor VAF data using mixtools
      
      mix_results <- normalmixEM(vafsByCn[[i]]$vaf, k=num_clusters, maxit=10000, maxrestarts=20);
      posteriors <- mix_results$posterior;
      clusters = NULL;
      for (n in 1:(dim(posteriors)[1])) {
        clusters[n]=as.numeric(which(posteriors[n,]==max(posteriors[n,]),arr.ind=T));
      }
      plot_clusters(vafsByCn[[i]]$vaf,clusters);

      ## print output file displaying data points clustered and their associated cluster
      if (!is.null(clusteredDataOutputFile)) {
        output = cbind(vafsByCn[[i]],clusters);
        names(output)[8] = "cluster"
        write.table(output, file=clusteredDataOutputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
      }
    }
  }
  axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.4);


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
drawScatterPlot <- function(data, highlightSexChrs, positionsToHighlight, colors, cn, maxDepth, highlightsHaveNames){

  ## define plot colors
  ptcolor = colors[cn]
  circlecolor = substr(colors[cn],1,7) #chop off the alpha value

  ## define the plot space by plotting offscreen points
  plot.default( x=-10000, y=1, log="y", type="p", pch=19, cex=0.4,
               col="#00000000", xlim=c(0,100), ylim=c(5,maxDepth*3),
               axes=FALSE, ann=FALSE, xaxs="i", yaxs="i");

  if(length(data[,1]) > 0){
    ## plot data points
    if(highlightSexChrs){
      ##plot autosomes
      data.autosomes = data[!(data$chr == "X" | data$chr == "Y"),]
      points(data.autosomes$vaf, data.autosomes$depth, type="p", pch=16, cex=0.75, col=ptcolor);
      ##plot sex chromsomes highlighted
      data.sex = data[(data$chr == "X" | data$chr == "Y"),]
      points(data.sex$vaf, data.sex$depth, type="p", pch=2, cex=0.75, col=ptcolor);
    } else {
      points(data$vaf, data$depth, type="p", pch=16, cex=0.75, col=ptcolor);
    }

    ##TODO - add a legend for point types - highlight vs sex chrs vs autosomes


    ## add highlighted of points selected for by user
    if(!(is.null(positionsToHighlight))){
      addpts = merge(data, positionsToHighlight, by.x=c(1,2), by.y = c(1,2))
      if(length(addpts[,1]) > 1){
        if(highlightsHaveNames){
          for(i in 1:length(addpts$vaf)){
            text(addpts$vaf[i],addpts$depth[i],labels=i,cex=0.5)
          }
        } else {
          points(x=addpts$vaf,y=addpts$depth,type="p",pch=7,cex=0.8,col="#555555FF");
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
##      ncol=ceiling(length(addpts[,1])/13)
##      legend(0, 750, legend=addpts$name, pch=as.character(1:(length(addpts[,1]))),
##             cex=0.5, bty="n", ncol=ncol)
      ypos=rev(seq(0,900,(900/13)))[1:13]
      ncol=ceiling(length(addpts[,1])/13)
      xpos=0;
      offset=1
      for(i in 1:ncol){        
        names = addpts[offset:(offset+12),]$name;
        names = as.character(names[!(is.na(names))])
        num = length(names)
        print(num)
        print(xpos)
        print(ypos[1:num])
        print(names)

        for(i in 1:num){
          text(xpos, ypos[i], paste(offset+i-1,". ",names[i],sep=""), cex=0.5, pos=4)
        }
        xpos=xpos+250;
        offset=offset+13;
      }
    }
  }
}


##--------------------------------------------------------------------
## function to plot clusters
plot_clusters <- function(data,clusters) {
point_colors <- rainbow(max(clusters),alpha=0.2);
  text_colors <- rainbow(max(clusters),alpha=1.0);

  for(k in 1:max(clusters)) {
    vector_size = rep(7.1,length(subset(data,clusters==k)));
    if (length(subset(data,clusters==k)) > 0) { #sometimes it doesn't assign any points to a cluster...don't understand why
      points(subset(data,clusters==k),vector_size,col=point_colors[k],cex=0.75,type="p",pch=16);
      ##text(x=mean(subset(data,clusters==k)),y=10,col=text_colors[k],labels=paste("C",k),cex=0.4);
      text(x=median(subset(data,clusters==k)),y=14.1,col=text_colors[k],labels=round(mean(subset(data,clusters==k)),1),cex=0.7);
      abline(v = (min(subset(data,clusters==k))-.4),col=text_colors[k],cex=0.5,lwd=0.4) # add vertical line at min of x
      abline(v = (max(subset(data,clusters==k))+.4),col=text_colors[k],cex=0.5,lwd=0.4) # add vertical line at max of x
      #print(min(subset(data,clusters==k)));
      #print(max(subset(data,clusters==k)));
    }
  }
}


