##--------------------------------------------------------------------
## hand off the data to the appropriate clustering algorithm
##
## Each clustering algorithm function should take a NxM matrix of vafs as input
## (where M = number of samples and N = number of points to cluster)
## and return a list containing the following items
##   cluster.assignments = a vector of length N containing a numeric cluster assignment
##   cluster.means = a matrix of size M x number_of_clusters containing the mean values
##                   of a cluster's vafs in each sample
##   cluster.upper = same as cluster.means, but containing upper confidence bounds instead of mean
##   cluster.lower = same as cluster.means, but containing lower confidence bounds instead of mean
##   fit.x = a matrix of size M x P, where P is an arbitrary number of points between 0 and 1
##           where the model fit was sampled
##   fit.y = a matrix of the same size as fit.x, containing the corresponding Y value for each X
##           Y values should be scaled between 0 and 1

clusterVafs <- function(vafMatrix, method="bmm", purities=100, params=NULL, samples=1){
  ##check for suitable method
  if(method == "bmm"){
   return(clusterWithBmm(vafMatrix, samples=samples))
  ## }  else if(method != "mixtoolsBinomial"){
  ##   return(clusterWithMixtools(vafs, "Binomial", purity, params));
  ## } else if (method != "mixtoolsNormal"){
  ##   return(clusterWithMixtools(vafs, "Normal", purity, params));
  } else {
    print("Error: please choose a supported clustering method\n[bmm|kmeans]");
    return(0);
  }
}



##--------------------------------------------------------------------------
## Do clustering with bmm (binomial mixture model)
##
clusterWithBmm <- function(vafs, initialClusters=10, samples=1){
  library(bmm)

  #replace any values of zero with a very small number to prevent errors
  delta <- .Machine$double.eps
  vafs[which(vafs==0)] = delta

  ## Initialize the hyperparameters of the Beta mixture model (bmm).
  hyperparams <- init.bmm.hyperparameters(vafs, initialClusters)

  ## Initialize the parameters of the bmm.
  params <- init.bmm.parameters(vafs, initialClusters, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0)

  ## Perform the clustering.
  ## Start with the provided number of clusters, but prune any with low probability
  bmm.res <- bmm(vafs, initialClusters, params$r, params$mu, params$alpha, params$nu,
                 params$beta, params$c, hyperparams$mu0, hyperparams$alpha0,
                 hyperparams$nu0, hyperparams$beta0, hyperparams$c0,
                 convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
  if(bmm.res$retVal != 0) {
    cat("WARNING: bmm failed to converge. No clusters assigned\n")
    return(NULL);
  }

  ##get the assignment of each point to a cluster
  probs = exp(1)^bmm.res$ln.rho
  numPoints = length(probs[,1])
  numClusters = length(probs[1,])
  clusters <- rep(NA,numPoints)
  for(n in 1:numPoints) {
    max.cluster <- 0
    max.assignment <- -1
    for(k in 1:numClusters) {
      if ( probs[n,k] > max.assignment ) {
        max.assignment <- probs[n,k]
        max.cluster <- k
      }
    }
    clusters[n] <- max.cluster
  }

  ## find confidence intervals around the means of the clusters
  intervals = bmm.narrowest.mean.interval.about.centers(bmm.res$mu, bmm.res$alpha, bmm.res$nu, bmm.res$beta, 0.68)
  means = intervals$centers
  lower = intervals$lb
  upper = intervals$ub

  
  ## Generate (x,y) values of the posterior predictive density
  n <- 1000

  ## Don't evaluate at x=0 or x=1, which will blow up
  x <- seq(0, 1, 1/n)[2:n]

  ## create a num_dimensions x n matrix of y values
  n=n-1;
  y <- rep.int(0, n)
  y = t(matrix(rep(y,dim(vafs)[2]),ncol=dim(vafs)[2]))

  ##for each dimension
  for(dim in 1:dim(vafs)[2]){
    ym <- matrix(data=0, nrow=numClusters, ncol=n)
    num.iterations <- 100
    for (k in 1:numClusters) {
      for (i in 1:n) {
        ## Evaluate posterior probability at x.
        ym[k,i] <- bmm.component.posterior.predictive.density(x[i], bmm.res$mu[dim,k], bmm.res$alpha[dim,k], bmm.res$nu[dim,k], bmm.res$beta[dim,k], bmm.res$E.pi[k], num.samples = num.iterations)
        y[dim,i] <- y[dim,i] + ym[k,i]
      }
    }
    ##scale yvals between 0 and 1
    y[dim,] = y[dim,]/max(y[dim,])
  }

  ##scale xvals between 1 and 100
  x = x*100
    
  #return a list of info
  return(list(
              cluster.assignments = clusters,
              cluster.means = means,
              cluster.upper = upper,
              cluster.lower = lower,
              fit.x = x,
              fit.y = y))
}


## ##--------------------------------------------------------------------------
## ## Do clustering with mixtools/mixdist
## ##
## clusterWithMixtools <- function(vafs,componentDistribution="Binomial"){
##   library(mixtools)
##   library(mixdist)


##   ##TODO set default params here
##   maximumClustersToTest=10;

##   ## ##overwrite params if any passed in
##   ## if(!(is.null(params))){
##   ##    
##   ## }
  
##   num_clusters = 0
##   print("Performing 'mixdist' clustering...");
  
##   ##function to process mixdist results
##   process_percents <- function(percents,chisq,pval,componentDistribution) {
##     minchisq = NULL;
##     best_fit_component_assessment = 0;
##     minpval = NULL;
##     true_cluster_count = NULL;

##     for (i in 1:maximumClustersToTest) {
##       if (percents[i]!="Error" && !is.nan(pval[i]) && pval[i] < 0.05) {
##         if (is.null(minchisq) || chisq[i] < minchisq) {
##           minchisq = chisq[i];
##           minpval = pval[i];
##           best_fit_component_assessment = i;
##         }
##       }
##       true_cluster_count[i] = 0;
##       percentage_per_cluster = as.numeric(percents[[i]]);
##       true_clusters = percentage_per_cluster > 0.02;
##       for (j in 1:i) {
##         if (isTRUE(true_clusters[j])) {
##           true_cluster_count[i] = true_cluster_count[i] + 1;
##         }
##       }
##     }
##     print(paste("chosen_component_assessment = ",best_fit_component_assessment,", true_component_count = ",true_cluster_count[best_fit_component_assessment],", assessment_chisq_value = ",minchisq,", assessment_p-value = ",minpval,sep=""));

##     ## If no clusters exceeded the pval of 0.05,it returns numeric(0),
##     ## which causes problems later on solution is to return just 0
##     ## clusters and warn user
##     if(best_fit_component_assessment > 0){
##       return(true_cluster_count[best_fit_component_assessment]);
##     } else {
##       print("WARNING: unable to estimate number of clusters - significance of 0.05 not reached")
##       return(0)
##     }
##   }
  
##   ## run mixdist
##   data = vafs$vaf/100
##   grouped_data = NULL;

##   ## if ceiling of max(data) is odd, then add 1 and group data. else, group data using the even ceiling
##   if (ceiling(max(data))%%2) {
##     grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data))+1,2)));
##   } else {
##     grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data)),2)));
##   }

##   ## for each component count, get mixdist fit estimates for chosen distribution
##   percents=NULL; chisq=NULL; pval=NULL;
##   for (i in 1:maximumClustersToTest) {
##     data_params = NULL;
##     test = NULL;

##     if (componentDistribution == "Normal") {
##       data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data),i));
##       test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="norm"), silent=TRUE);
##     }
##     if (componentDistribution == "Binomial") {
##       data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data)/2,i));
##       test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="binom", constr=mixconstr(consigma="BINOM",size=rep(round(length(data)),i))), silent=TRUE);
##     }

##     if (class(test) == 'try-error') {
##       percents[i] = "Error";
##       print(paste("Common mixdist error when looking for ",i," components.",sep=""));
##     }
##     else {
##       percents[i]=list(test$parameters$pi)
##       chisq[i]=test$chisq;
##       pval[i] = test$P;

##       if(testing){
##         filename = paste("plot_component_",i,".pdf",sep="");
##         print(paste("Testing output: plotting ",filename));
##         pdf(filename)
##         ##dev.new();
##         plot(test);
##         ##dev.copy(pdf,filename);
##         d = dev.off();
##       }
##     }
##   }
##   num_clusters = process_percents(percents,chisq,pval,componentDistribution);


##   ##now that we know the number of clusters, assign points to each using mixtools
##   print("assigning points to clusters using mixtools");
##   mix_results <- normalmixEM(vafs, k=num_clusters, maxit=10000, maxrestarts=20);
##   posteriors <- mix_results$posterior;
##   clusters = NULL;
##   for (n in 1:(dim(posteriors)[1])) {
##     clusters[n]=as.numeric(which(posteriors[n,]==max(posteriors[n,]),arr.ind=T));
##   }
##   plot_clusters(vafsByCn[[i]]$vaf,clusters);
  

##   ## print output file displaying data points clustered and their associated cluster
##   if (!is.null(clusteredDataOutputFile)) {
##     output = cbind(vafsByCn[[i]],clusters);
##     names(output)[8] = "cluster"
##     write.table(output, file=clusteredDataOutputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
##   }  
  

## }
