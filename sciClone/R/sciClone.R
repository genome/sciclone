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
sciClone <- function(vafs, outputPrefix, copyNumberCalls=NULL, regionsToExclude = NULL,
                     sampleNames, minimumDepth=100, clusterMethod="bmm",clusterParams=NULL,
                     minimumLabelledPeakHeight=0.001, onlyLabelHighestPeak=FALSE,
                     overlayClusters=FALSE, overlayIndividualModels=FALSE,
                     show1DHistogram=FALSE, plotOnlyCN2=FALSE,
                     showCopyNumberScatterPlots=TRUE, positionsToHighlight=NULL,
                     purities=NULL, highlightSexChrs=TRUE, cnCallsAreLog2=FALSE,
                     useSexChrs=TRUE, highlightsHaveNames=FALSE, doClustering=TRUE,
                     showTitle=TRUE, plotIntermediateResults=0, verbose=TRUE){

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
    vafs[[i]] = cleanAndAddCN(vafs[[i]], copyNumberCalls[[i]], i, cnCallsAreLog2, regionsToExclude, useSexChrs, minimumDepth)
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

  ## merge the data frames to get a df with vafs and readcounts for each variant in each sample
  vafs.merged = vafs[[1]]
  if(dimensions > 1){
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

  nonvafcols <- (1:length(names(vafs.merged)))[!((1:length(names(vafs.merged))) %in% vafcols)]

  vafs.merged.orig <- vafs.merged
  vafs.merged.cn2.orig <- vafs.merged.cn2
  
  #convert to a matrix to feed into clustering
  vafs.matrix = as.matrix(vafs.merged.cn2[,vafcols])
  #convert vafs to be between 0 and 1
  vafs.matrix = vafs.matrix/100

  clust=NULL
  if(doClustering){
    if(verbose){print("clustering...")}
    clust=clusterVafs(vafs.merged.cn2, vafs.matrix, clusterMethod, purities, clusterParams, samples=length(purities), plotIntermediateResults=plotIntermediateResults)
    if(verbose){print("finished clustering full-dimensional data...");}
  }

  numClusters=0
  if(!(is.null(clust))){
    numClusters = max(clust$cluster.assignments,na.rm=T)
    #append (hard and fuzzy) cluster assignments

    vafs.merged.cn2 = cbind(vafs.merged.cn2,cluster=clust$cluster.assignments)
    vafs.merged.cn2 = cbind(vafs.merged.cn2,cluster.prob=clust$cluster.probabilities)
    vafs.merged = merge(vafs.merged,vafs.merged.cn2, by.x=c(1:length(vafs.merged)), by.y=c(1:length(vafs.merged)),all.x=TRUE)
    #sort by chr, st
    vafs.merged = vafs.merged[order(vafs.merged[,1], vafs.merged[,2]),]
    print(paste("found",numClusters,"clusters"))
  }

  # Show a 1D projection of the data along the axes for multidimensional data.
  showMarginalData <- TRUE
  #showMarginalData <- FALSE
  if(dimensions == 1) { showMarginalData <- FALSE }
  # Cluster the 1D data shown along the margins (as opposed to showing the
  # 1D projection of the multidimensional clustering result)
  doClusteringAlongMargins <- TRUE
  if(doClustering == FALSE) { doClusteringAlongMargins <- FALSE }
  if(dimensions == 1) { doClusteringAlongMargins <- FALSE }
  if(showMarginalData == FALSE) { doClusteringAlongMargins <- FALSE }

  # Perform 1D clustering of each dimension independently.
  marginalClust = list()
  vafs.1d = list()
  if(doClusteringAlongMargins == TRUE){
    for(i in 1:dimensions){
      marginalClust[[i]]=clusterVafs(NULL, vafs.matrix[,i,drop=FALSE], clusterMethod, purities[i], clusterParams, FALSE)
      if(verbose){print(paste("finished clustering", sampleNames[i], "..."))}
      numClusters = max(marginalClust[[i]]$cluster.assignments,na.rm=T)
      print(paste("found",numClusters,"clusters"))

      vafs.1d.merged.cn2 = cbind(vafs.merged.cn2.orig,cluster=marginalClust[[i]]$cluster.assignments)
      vafs.1d.merged = merge(vafs.merged.orig,vafs.1d.merged.cn2, by.x=c(1:length(vafs.merged.orig)), by.y=c(1:length(vafs.merged.orig)),all.x=TRUE)
      #sort by chr, st
      vafs.1d.merged = vafs.1d.merged[order(vafs.1d.merged[,1,drop=FALSE], vafs.1d.merged[,2,drop=FALSE]),]
      vafs.1d[[i]] = vafs.1d.merged
    }
  }
  
  ##--------------------------------------------------------
  ##output and plots
  clusteredDataOutputFile=paste(outputPrefix,".clusters",sep="");
  if(!(is.null(clust)) & !(is.null(clusteredDataOutputFile))){
    ##TODO - order these first    
    write.table(vafs.merged, file=clusteredDataOutputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
  }
  
  plot1d(vafs.merged, outputPrefix, densityData, sampleNames, dimensions, plotOnlyCN2,
         showCopyNumberScatterPlots, clust, highlightSexChrs, positionsToHighlight,
         highlightsHaveNames, overlayClusters, overlayIndividualModels, show1DHistogram,
         onlyLabelHighestPeak, minimumLabelledPeakHeight, showTitle);

  if((showMarginalData == TRUE) & (doClusteringAlongMargins == TRUE)) {
    plot2dWithMargins(vafs.1d, vafs.merged, outputPrefix, densityData, sampleNames, dimensions, plotOnlyCN2, 
                      marginalClust, clust, highlightSexChrs, positionsToHighlight, highlightsHaveNames,
                      overlayClusters, onlyLabelHighestPeak, minimumLabelledPeakHeight)
  }
  
  if(dimensions>1) {
    plot2d(vafs.merged, outputPrefix, sampleNames, dimensions, positionsToHighlight, highlightsHaveNames, overlayClusters)
  }
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

  ##round these values to absolute calls
  for(n in 0:4){
    pos = which(vafs$cn >= (n-0.5) & vafs$cn < n+0.5)
    if(length(pos) > 0){
      vafs[pos,]$cn = n
    }
  }
    
  if(length(which(is.na(vafs$cn))) > 0){
    print("Not all variants fall within a provided copy number region. The copy number of these variants is assumed to be 2.")
    vafs[which(is.na(vafs$cn)),]$cn = 2
  }

  return(vafs)
}


##--------------------------------------------------------------------
## intersect the variants with the regionsToExclude to remove them
##
excludeRegions <- function(vafs,regionsToExclude){
    library(IRanges);
    ##for each chromosome, find variants falling inside regions to be excluded
    for(chr in names(table(regionsToExclude[,1]))){
        vars = IRanges(start=vafs[vafs$chr==chr,]$st, end=vafs[vafs$chr==chr,]$st);
        excludedRegions = IRanges(start=regionsToExclude[regionsToExclude[,1]==chr,2], end=regionsToExclude[regionsToExclude[,1]==chr,3]);
        if( (length(vars) == 0) | (length(excludedRegions) == 0)){
            next;
        }
        vars_to_exclude = as.matrix(findOverlaps(vars,excludedRegions));
        ##if there are variants that fell inside the exclude regions, find them and remove them from vafs
        if(length(vars_to_exclude) > 0){
            for(i in vars_to_exclude[,1]){
                vpos = which(vafs$chr==chr & vafs$st==start(vars[i,]));
                if(identical(vpos,integer(0))) { next; }
                ##remove the variant from vafs
                vafs = vafs[-vpos,];
            }
        }
    }

    return(vafs)
} # end excludeRegions


##---------------------------------------------------------------------
## clean up vaf data, add cn
##
cleanAndAddCN <- function(vafs, cn, num, cnCallsAreLog2, regionsToExclude, useSexChrs, minimumDepth){
    names(vafs) = c("chr","st","ref","var","vaf")

    ##remove MT values
    vafs = vafs[!(vafs$chr == "M" | vafs$chr == "MT"),]
    if(length(vafs[,1]) == 0){return(vafs)}
    
    ##remove NA sites
    vafs = vafs[!(is.na(vafs$vaf)),]
    if(length(vafs[,1]) == 0){return(vafs)}
    
    ##remove sites in excludedRegions
    vafs = excludeRegions(vafs,regionsToExclude);
    if(length(vafs[,1]) == 0){return(vafs)}
    
    ##add depth
    vafs = vafs[vafs$vaf > 0,]
    vafs$depth = round(vafs$var/(vafs$vaf/100))

    ##add cn calls
    if(is.null(cn)){
        ##assume all sites are 2x if no cn info
        vafs$cn = 2;
    } else {
        if(cnCallsAreLog2){
            cn[,4] = (2^(cn[,4]))*2
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


