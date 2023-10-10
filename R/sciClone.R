##---------------------------------------------------------------
## clean up and cluster variants based on allele frequency
##
sciClone <- function(vafs, copyNumberCalls=NULL, regionsToExclude=NULL,
                     sampleNames, minimumDepth=100, clusterMethod="bmm",
                     clusterParams="no.apply.overlapping.std.dev.condition",
                     cnCallsAreLog2=FALSE,
                     useSexChrs=TRUE, doClustering=TRUE, verbose=TRUE,
                     copyNumberMargins=0.25, maximumClusters=10, annotation=NULL,
                     doClusteringAlongMargins=TRUE, plotIntermediateResults=0){

  if(verbose){print("checking input data...")}

  #how many samples do we have?
  dimensions=NULL;
  if(is.data.frame(vafs)){
    dimensions = 1;
    vafs=list(vafs)
    copyNumberCalls=list(copyNumberCalls)
  } else if(is.list(vafs)){
    dimensions = length(vafs);
  } else {
    stop("input param vafs must be either a data frame (for 1-sample clustering), or a list of data frames (for multi-sample clustering)")
  }

  if(missing(sampleNames)){
    stop("sampleNames is a required parameter")
  }
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

  #roughly implemented, but not reliable yet, so always turned off at present
  doPurityScaling=FALSE
  purities=NULL;


  if(!(is.null(purities))){
    if(length(purities) != dimensions){
      stop(paste("the number of input purities calls(",length(purities),") does not equal the number of input samples (",dimensions,")\nEither provide a purity for each sample, or set purities to NULL and it will be estimated for you",sep=""))
    }
  }

  if(!is.null(regionsToExclude)){
    if(is.data.frame(regionsToExclude)){
      regionsToExclude = list(regionsToExclude);
    }
  }




  ##-----------------------------------------
  densityData=NULL
  ##clean up data, get kernel density, estimate purity
  for(i in 1:dimensions){
    vafs[[i]] = cleanAndAddCN(vafs[[i]], copyNumberCalls[[i]], i, cnCallsAreLog2, regionsToExclude, useSexChrs, minimumDepth, copyNumberMargins)

    ##calculate the densities and peaks for variants of each copy number
    if(is.null(densityData)){
      densityData = list(getDensity(vafs[[i]], copyNumberMargins, minimumDepth))
    } else {
      densityData= c(densityData, list(getDensity(vafs[[i]], copyNumberMargins, minimumDepth)))
    }

    ## This stop should probably be replaced so that plotting can take place
    ## for tumors with ployplody, even if we can't cluster
    ## (maybe we should even cluster with 3x regions, etc - put it on the todo list)
    if(is.null(densityData[[i]]$densities[[2]])){
      cat(paste("can't do clustering - no copy number 2 regions to operate on in sample",i,"\n"));
      if(doClustering==TRUE) { return(NULL) }
    }
  }



  ##-----------------------------------------------
  ## munge the data

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
  cleancncols = grep("^cleancn",names(vafs.merged))

  ##change NA values introduced by merge to ref/var/vaf/depth of zero, cn of 2
  for(i in c(vafcols,refcols,varcols,depthcols)){
    vafs.merged[is.na(vafs.merged[,i]),i] = 0;
  }

  #add sample names to make output pretty
  names(vafs.merged)[refcols] = paste(sampleNames,".ref",sep="")
  names(vafs.merged)[varcols] = paste(sampleNames,".var",sep="")
  names(vafs.merged)[vafcols] = paste(sampleNames,".vaf",sep="")
  names(vafs.merged)[depthcols] = paste(sampleNames,".depth",sep="")
  names(vafs.merged)[cncols] = paste(sampleNames,".cn",sep="")
  names(vafs.merged)[cleancncols] = paste(sampleNames,".cleancn",sep="")

  #determine whether there is adequate depth at each site in all samples
  adequateDepth =  rep(1,length(vafs.merged[,1]))
  for(i in 1:length(vafs.merged[,1])){
    if(length(which(vafs.merged[i,depthcols] >= minimumDepth)) < length(depthcols)){
      adequateDepth[i] = 0
    }
  }
  vafs.merged$adequateDepth=adequateDepth

  #determine whether all samples are copy number neutral at each site
  cnNeutral = rep(1,length(vafs.merged[,1]))
  for(i in 1:length(vafs.merged[,1])){
    for(j in cleancncols){
      if(is.na(vafs.merged[i,j])){
        cnNeutral[i] = 0
      } else if(vafs.merged[i,j] != 2){
        cnNeutral[i] = 0
      }
    }
  }

  ## remove any lines where all clean CN columns are not 2 and depth is adequate
  ## we only cluster based on sites that are CN neutral in all samples
  vafs.merged.cn2 = vafs.merged[(as.logical(cnNeutral) & as.logical(adequateDepth)),]
  if(length(vafs.merged.cn2[,1]) < 1){
    print("ERROR: no sites are copy number neutral and have adequate depth in all samples")
    return(NULL);
  }

  cat(paste(length(vafs.merged.cn2[,1]),"sites (of", length(vafs.merged[,1]), "original sites) are copy number neutral and have adequate depth in all samples\n"))
  cat(paste(nrow(vafs.merged[!as.logical(cnNeutral),])), "sites (of", nrow(vafs.merged), "original sites) were removed because of copy-number alterations\n")
  cat(paste(nrow(vafs.merged[!as.logical(adequateDepth),])), "sites (of", nrow(vafs.merged), "original sites) were removed because of inadequate depth\n")
  cat(paste(nrow(vafs.merged[!as.logical(cnNeutral) | !as.logical(adequateDepth),])), "sites (of", nrow(vafs.merged), "original sites) were removed because of copy-number alterations or inadequate depth\n")

  nonvafcols <- (1:length(names(vafs.merged)))[!((1:length(names(vafs.merged))) %in% vafcols)]

  vafs.merged.orig <- vafs.merged
  vafs.merged.cn2.orig <- vafs.merged.cn2

  #convert to a matrix to feed into clustering
  vafs.matrix = as.matrix(vafs.merged.cn2[,vafcols])
  vars.matrix = as.matrix(vafs.merged.cn2[,varcols])
  refs.matrix = as.matrix(vafs.merged.cn2[,refcols])
  #convert vafs to be between 0 and 1
  vafs.matrix = vafs.matrix/100

  
  ## purity correction - always turned off until it can made more reliable.
  if(is.null(purities)){
    if(doPurityScaling){
      purities=c()
      print("estimating purity")

      ##do multi-d clustering of the data, use the higest point in the max cluster to estimate purity (outliers are removed in clustering)
      clust=clusterVafs(vafs.merged.cn2, vafs.matrix, vars.matrix, refs.matrix, maximumClusters, clusterMethod, clusterParams,
        samples=length(purities), plotIntermediateResults=0, verbose=0)
      if(is.null(clust[[1]])) {
        print("WARNING: couldn't estimate and correct for purity, will cluster using input vafs")
      } else {
        for(i in 1:dimensions){
          ##use the mean of the max cluster
          purities[i] = round(max(sort(t(clust[["cluster.means"]])[,i]))*2*100,2)

          #use the highest point that was assigned to a cluster (not an outlier)
          vafs.tmp = cbind(vafs.merged.cn2,cluster=clust$cluster.assignments)
          vafcols = grep("vaf$",names(vafs.tmp))
          #purities[i] = max(vafs.tmp[vafs.tmp$cluster > 0,vafcols[i]])*2
          print(paste("purity of sample",sampleNames[i],"is estimated to be",purities[[i]]))
          ##do the adjustment to the appropriate columns
          vafcols = grep("vaf$",names(vafs.merged))
          vafs.merged[,vafcols[i]] = (vafs.merged[,vafcols[i]]/purities[i])*100
          vafcols = grep("vaf$",names(vafs.merged.cn2))
          vafs.merged.cn2[,vafcols[i]] = (vafs.merged.cn2[,vafcols[i]]/purities[i])*100
          vafs.matrix[,i] = (vafs.matrix[,i]/purities[i])*100
        }
      }

    } else {
      purities = rep(100,dimensions)
    }
  } else { ##purities provided
    if(doPurityScaling){
      ##use input purities for adjustment
      for(i in 1:dimensions){
        vafs.merged[,vafcols[i]] = (vafs.merged[,vafcols[i]]/purities[i])*100
        vafs.merged.cn2[,vafcols[i]] = (vafs.merged.cn2[,vafcols[i]]/purities[i])*100
        vafs.matrix[,i] = (vafs.matrix[,i]/purities[i])*100
      }
    }
  }

  ##---------------------------------------------------
  ##do the clustering
  marginalClust = list()
  vafs.1d = list()
  if(dimensions == 1) { doClusteringAlongMargins <- FALSE }
  if(doClustering == FALSE) { doClusteringAlongMargins <- FALSE }

  ## Perform 1D clustering of each dimension independently.
  if(doClusteringAlongMargins == TRUE){
    print("clustering each dimension independently")
    for(i in 1:dimensions){
      marginalClust[[i]]=clusterVafs(NULL, vafs.matrix[,i,drop=FALSE], vars.matrix[,i,drop=FALSE], refs.matrix[,i,drop=FALSE],
                     maximumClusters, clusterMethod, clusterParams, FALSE)
      if(verbose){print(paste("finished 1d clustering", sampleNames[i], "..."))}
      numClusters = max(marginalClust[[i]]$cluster.assignments,na.rm=T)
      print(paste("found",numClusters,"clusters using", clusterMethod, "in dimension",sampleNames[i]))
      print(marginalClust[[i]]$cluster.means)

      vafs.1d.merged.cn2 = cbind(vafs.merged.cn2.orig,cluster=marginalClust[[i]]$cluster.assignments)
      vafs.1d.merged = merge(vafs.merged.orig,vafs.1d.merged.cn2, by.x=c(1:length(vafs.merged.orig)),
        by.y=c(1:length(vafs.merged.orig)),all.x=TRUE)
      ##sort by chr, st
      # vafs.1d.merged = vafs.1d.merged[order(vafs.1d.merged[,1,drop=FALSE], vafs.1d.merged[,2,drop=FALSE]),]
      # hot fix using arrange from dplyr package (axel kunstner 2023-10-10)
      vafs.1d.merged = vafs.1d.merged |> dplyr::arrange(chr, st)
      vafs.1d[[i]] = vafs.1d.merged
    }
  }

  ## now do clustering using all dimensions
  clust=list(NULL)
  if(doClustering){
    if(verbose){print("clustering...")}
    clust=clusterVafs(vafs.merged.cn2, vafs.matrix, vars.matrix, refs.matrix, maximumClusters, clusterMethod, clusterParams,
      samples=length(purities), plotIntermediateResults=plotIntermediateResults, verbose=0)
    if(is.null(clust[[1]])) {
      print("Warning: no clusters, returning NULL")
      return(NULL)
    }
    if(verbose){print("finished clustering full-dimensional data...");}

    ##reorder the clusters to the largest vaf is first
    clust=reorderClust(clust)

  }


  numClusters=0
  if(!(is.null(clust[[1]]))){
    numClusters = max(clust$cluster.assignments,na.rm=T)
    ##append (hard and fuzzy) cluster assignments
    vafs.merged.cn2 = cbind(vafs.merged.cn2,cluster=clust$cluster.assignments)
    vafs.merged.cn2 = cbind(vafs.merged.cn2,cluster.prob=clust$cluster.probabilities)
    vafs.merged = merge(vafs.merged,vafs.merged.cn2, by.x=c(1:length(vafs.merged)), by.y=c(1:length(vafs.merged)),all.x=TRUE)
    ##sort by chr, st
    vafs.merged = vafs.merged[order(vafs.merged[,1], vafs.merged[,2]),]
    print(paste("found",numClusters,"clusters using", clusterMethod, "in full dimensional data"))
  }

  if(!is.null(annotation)) {
    vafs.merged = merge(vafs.merged, annotation, by.x=c(1,2), by.y=c(1,2), all.x=TRUE, all.y=FALSE)
  }
  return(new("scObject", clust=clust, densities=densityData, dimensions=dimensions,
             marginalClust=marginalClust, sampleNames=sampleNames, vafs.1d=vafs.1d,
             vafs.merged=vafs.merged, purities=purities))
}


##------------------------------------------------------------------------------------------
##  write out a table of all vafs, cns, and cluster assignments
##
writeClusterTable <- function(sco, outputFile){
    write.table(sco@vafs.merged, file=outputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
}

##------------------------------------------------------------------------------------------
##  write out a table of all cluster centers and their SEMs
##
writeClusterSummaryTable <- function(sco, outputFile){
    out <- paste(outputFile, ".means", sep="")
    a = t(sco@clust[["cluster.means"]])
    #num.clusters <- max(sco@clust$cluster.assignments)
    num.clusters <- nrow(a)
    colnames(a) = c(sco@sampleNames)
    rownames(a) = paste("cluster",1:num.clusters,sep="")
    write.table(a,file=out,row.names=TRUE,col.names=NA,sep="\t",quote=F)

    out <- paste(outputFile, ".lower", sep="")
    a = t(sco@clust[["cluster.lower"]])
    colnames(a) = c(sco@sampleNames)
    rownames(a) = paste("cluster",1:num.clusters,sep="")
    write.table(a,file=out,row.names=TRUE,col.names=NA,sep="\t",quote=F)

    out <- paste(outputFile, ".upper", sep="")
    a = t(sco@clust[["cluster.upper"]])
    colnames(a) = c(sco@sampleNames)
    rownames(a) = paste("cluster",1:num.clusters,sep="")
    write.table(a,file=out,row.names=TRUE,col.names=NA,sep="\t",quote=F)
}

##------------------------------------------------------------------------------------------
##  return an N x N connectivity matrix C where c_ij = 1 iff items i and j
##  are co-clustered, for all items i,j = 1 to N.
##
getConnectivityMatrix <- function(sco){
  vaf.table <- sco@vafs.merged
  N <- dim(vaf.table)[1]
  C <- matrix(data=0, ncol=N, nrow=N)
  clusters <- vaf.table$cluster
  for(i in 1:N) {
    for(j in 1:N) {
      if(!is.na(clusters[i]) & !is.na(clusters[j]) & (clusters[i] == clusters[j])) { C[i,j] <- 1 }
    }
  }
  return(C)
}

## -----------------------------------------------------
## Calculate the pvalue of a vaf v being in cluster k
calculate.pvalue <- function(vaf, k, sco) {
  cluster.method <- sco@clust$cluster.method
  if(cluster.method == "bmm") {
    mu <- sco@clust$mu
    alpha <- sco@clust$alpha
    nu <- sco@clust$nu
    beta <- sco@clust$beta
    pi <- sco@clust$pi
    pvalue <- bmm.calculate.pvalue(vaf, k, mu, alpha, nu, beta, pi)
  } else {
    stop(paste("calculate.pvalue not implemented for", cluster.method))
  }
  return(pvalue)
}

##--------------------------------------------------------------------
## intersect the variants with CN calls to classify them
##
addCnToVafs <- function(vafs, cncalls, copyNumberMargins){
  suppressPackageStartupMessages(library(IRanges))
  vafs$cn = NA
  vafs$cleancn = NA
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
        if(length(matches[which(matches[,1]==i),2])>1){
          print("ERROR: input site matches two copy number regions (this should be impossible):")
          print(cnsegs[matches[which(matches[,1]==i),2]])
          print("Is your CN file in one-based format?")
          stop("unable to continue")
        }
        cpos = which(cncalls[,1] == chr & cncalls[,2]==start(cnsegs[matches[which(matches[,1]==i),2]]))
        ##set the value
        vafs[vpos,]$cn = cncalls[cpos,4]
      }
    }
  }

  ##round these values to absolute calls
  for(n in 0:4){
    pos = which(vafs$cn >= (n-copyNumberMargins) & vafs$cn < (n+copyNumberMargins))
    if(length(pos) > 0){
      vafs[pos,]$cleancn = n
    }
  }

  indices <- which(is.na(vafs$cn))
  if(length(indices) > 0){
    print("Not all variants fall within a provided copy number region. The copy number of these variants is assumed to be 2.")
    vafs[indices,]$cn = 2
    vafs[indices,]$cleancn = 2
  }

  return(vafs)
}

##--------------------------------------------------------------------
## intersect the variants with the regionsToExclude to remove them
##
excludeRegions <- function(vafs,regionsToExclude){

  #read in all the excluded regions and combine them into one file;
  regs = NULL
  for(i in 1:length(regionsToExclude)){
    a = regionsToExclude[[i]][,1:3];
    if(!is.null(regs)){
      a = rbind(regs,a)
    } else{
      regs = a
    }
  }
  #sort the list
  regs = regs[order(regs[,1], regs[,2]),]


  suppressPackageStartupMessages(library(IRanges))
  ##for each chromosome, find variants falling inside regions to be excluded
  for(chr in names(table(regs[,1]))){
    vars = IRanges(start=as.numeric(as.character(vafs[vafs$chr==chr,]$st)),
      end=as.numeric(as.character(vafs[vafs$chr==chr,]$st)));
    excludedRegions = IRanges(start=regs[regs[,1]==chr,2], end=regs[regs[,1]==chr,3]);
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
cleanAndAddCN <- function(vafs, cn, num, cnCallsAreLog2, regionsToExclude, useSexChrs, minimumDepth, copyNumberMargins){
    names(vafs) = c("chr","st","ref","var","vaf")

    ##remove MT values
    vafs = vafs[!(vafs$chr == "M" | vafs$chr == "MT"),]
    if(length(vafs[,1]) == 0){return(vafs)}

    ##remove NA sites
    vafs = vafs[!(is.na(vafs$vaf)),]
    if(length(vafs[,1]) == 0){return(vafs)}

    ##remove duplicate sites
    vafs = unique(vafs)

    ##make sure columns are numeric
    for(i in 2:5){
      if(!(is.numeric(vafs[,i]))){
        print(paste("ERROR: column",names(vafs)[i]," in sample",i,"is not numeric"))
        stop();
      }
    }

    ##remove sites in excludedRegions
    if(!is.null(regionsToExclude)){
      vafs = excludeRegions(vafs,regionsToExclude);
      if(length(vafs[,1]) == 0){return(vafs)}
    }

    ##add cn calls
    if(is.null(cn)){
        ##assume all sites are 2x if no cn info
        vafs$cn = 2;
        vafs$cleancn = 2;
    } else {
        if(cnCallsAreLog2){
            cn[,4] = (2^(cn[,4]))*2
        }
        vafs = addCnToVafs(vafs, cn, copyNumberMargins)
    }

    ##add depth
    vafs = vafs[vafs$vaf > 0 | ( vafs$var + vafs$ref ) > 0,]
    vafs$depth = mapply(function(var, ref, vaf) ifelse(vaf == 0, var + ref, round(var/(vaf/100))), vafs$var, vafs$ref, vafs$vaf)

    ##remove sex chromosomes if specified
    if(!(useSexChrs)){
        vafs = vafs[vafs$chr != "X" & vafs$chr != "Y",]
    }

    return(vafs)
}



##--------------------------------------------------------------------------
## calculate a sample's purity from VAF peaks - experimental,
## doesn't work all that well at the moment
##
getPurity <- function(peakPos){
  purity = 0
    if(length(peakPos[[2]][peakPos[[2]] <= 60]) > 0){
        purity = max(peakPos[[2]][peakPos[[2]] <= 50])*2;
        if(length(peakPos[[2]][peakPos[[2]] > 50]) > 0){
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
getDensity <- function(vafs, copyNumberMargins, minimumDepth){
    ##data structure to hold info
    densities = vector("list",4)
    factors = vector("list",4)
    peakPos = vector("list",4)
    peakHeights = vector("list",4)
    maxDensity = 0
    maxDepth=0

    for(i in 1:4){
        ##grab only the variants in this copy number and with adequate depth
        v = vafs[(vafs$cn > (i-copyNumberMargins)) & (vafs$cn < (i+copyNumberMargins)) & (vafs$depth > minimumDepth),]
        if(length(v[,1]) > 0){

            ##need two points for density calc
            if(length(v[,1]) > 1){
                ##calculate the density
                densities[[i]] = density(v$vaf, from=0, to=100, na.rm=TRUE, adj=0.85)
                factors[[i]] = (length(v[,1])/length(vafs[,1]))*densities[[i]]$y
                ##find the peaks
                p = c(getPeaks(factors[[i]]),FALSE,FALSE)
                peakPos[[i]] = densities[[i]]$x[p]
                peakHeights[[i]] = factors[[i]][p]
                ##store the largest density for use in scaling the plots later
                if(max(factors[[i]]) > maxDensity){
                    maxDensity = max(factors[[i]])
                }
                #filter out peaks unless they hit 10% of the max height
                keep=which(peakHeights[[i]] > maxDensity*0.10)
                peakPos[[i]] = peakPos[[i]][keep]
                peakHeights[[i]] = peakPos[[i]][keep]
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
## reorder the clusters so the largest vaf is first
##
reorderClust <- function(clust){

  nums=unique(clust$cluster.assignments)[unique(clust$cluster.assignments)>0]

  #calc distance of cluster center from origin
  dist=c()
  for(i in nums){
    z = clust$cluster.means
    dist = c(dist, sqrt(sum(clust$cluster.means[,i]^2)))
  }

  df = data.frame(nums,dist=dist)
  df = cbind(df[rev(order(df$dist)),],new=1:length(nums))


  ass = list()
  prob=list()
  means=list()
  upper=list()
  lower=list()
  individual=list()

  for(i in 1:length(df[,1])){
    oldnum = df[i,1]
    ## store the cluster assignments
    ass[[i]] = which(clust$cluster.assignments==oldnum)
    ##store the cluster probabilities
    prob[[i]] = clust$cluster.probabilities[,oldnum]
    means[[i]] = clust$cluster.means[,oldnum]
    lower[[i]] = clust$cluster.lower[,oldnum]
    upper[[i]] = clust$cluster.upper[,oldnum]
    individual[[i]] = clust$individual.fits.y[[oldnum]]
    if(is.null(clust$individual.fits.y[[oldnum]])) {
      individual[[i]] = 0
    }
  }
  for(i in 1:length(df[,1])){
    clust$cluster.assignments[ass[[i]]] = i
    clust$cluster.probabilities[,i] = prob[[i]]
    clust$cluster.means[,i] = means[[i]]
    clust$cluster.lower[,i] = lower[[i]]
    clust$cluster.upper[,i] = upper[[i]]
    clust$individual.fits.y[[i]] = individual[[i]]
  }

  if(clust$cluster.method == "bmm") {
    clust <- reorderBetaClust(clust, df)
  } else if (clust$cluster.method == "gaussian.bmm") {
    clust <- reorderGaussianClust(clust, df)
  } else if (clust$cluster.method == "binomial.bmm") {
    clust <- reorderBinomialClust(clust, df)
  } else {
    stop(paste("Cluster reordering not implemented for method", clust$cluster.method))
  }
  return(clust)
}
