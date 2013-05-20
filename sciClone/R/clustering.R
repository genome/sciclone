##--------------------------------------------------------------------
## hand off the data to the appropriate clustering algorithm
##
## Each clustering algorithm function should take as input
## an NxM matrix of vafs (with M = # of samples and N = # of points to cluster),
## an NxM matrix of variant read counts,
## and an NxM matrix of reference read counts.
## e.g., a beta mixture model works on the vafs, while a binomial mixture
## model operates on the variant and reference counts directly.  These
## two specifications are equivalent--we should probably drop the vaf input
## and have the method compute those if it needs them.
## It should return a list containing the following items
##   cluster.assignments = a vector of length N containing a numeric cluster assignment
##   cluster.means = a matrix of size M x number_of_clusters containing the mean values
##                   of a cluster's vafs in each sample
##   cluster.upper = same as cluster.means, but containing upper confidence bounds instead of mean
##   cluster.lower = same as cluster.means, but containing lower confidence bounds instead of mean
##   fit.x = a vector of length P, where P is an arbitrary number of points between 0 and 1
##           where the model fit was sampled (in each of the M dimensions)
##   fit.y = a matrix of size MxP, containing the corresponding Y value for each X
##           Y values should be scaled between 0 and 1
##   individual.fits.y = a list of length number_of_clusters, holding the individual fits for each of the models, each represented by an M x P matrix (as for fit.y)

clusterVafs <- function(vafs.merged, vafMatrix, varMatrix, refMatrix, maximumClusters, method="bmm", params=NULL, samples=1, plotIntermediateResults = 0, verbose=0){
  ##check for suitable method
  if(method == "bmm"){
    if(!is.null(params)){
      ##handle input params to clustering method - only supports one for now...
      params=strsplit(params,", *",perl=TRUE)[[1]]
      if(grepl("overlap.threshold",params)){
        val = strsplit(params[grep("overlap.threshold",params)] ," *= *",perl=TRUE)[[1]][2]
        return(clusterWithBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0, overlap.threshold=val,initialClusters=maximumClusters))
      }
    }
    return(clusterWithBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0,initialClusters=maximumClusters))
  } else if(method == "binomial.bmm"){
    if(!is.null(params)){
      ##handle input params to clustering method - only supports one for now...
      params=strsplit(params,", *",perl=TRUE)[[1]]
      if(grepl("overlap.threshold",params)){
        val = strsplit(params[grep("overlap.threshold",params)] ," *= *",perl=TRUE)[[1]][2]
        return(clusterWithBinomialBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0, overlap.threshold=val,initialClusters=maximumClusters))
      }
    }
    return(clusterWithBinomialBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0,initialClusters=maximumClusters))
  } else if(method == "gaussian.bmm"){
    if(!is.null(params)){
      ##handle input params to clustering method - only supports one for now...
      params=strsplit(params,", *",perl=TRUE)[[1]]
      if(grepl("overlap.threshold",params)){
        val = strsplit(params[grep("overlap.threshold",params)] ," *= *",perl=TRUE)[[1]][2]
        return(clusterWithGaussianBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0, overlap.threshold=val,initialClusters=maximumClusters))
      }
    }
    return(clusterWithGaussianBmm(vafs.merged, vafMatrix, varMatrix, refMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0,initialClusters=maximumClusters))
  } else {
    print("Error: please choose a supported clustering method\n[bmm|gaussian.bmm|binomial.bmm]");
    return(0);
  }
}


##--------------------------------------------------------------------------
## Go from fuzzy probabilities to hard cluster assignments
##
hardClusterAssignments <- function(numPoints,numClusters,probabilities) {
    # Any point that has been removed from the data set will have a
    # probability of NA and will be given an assignment of 0,
    # indicating that it is an outlier.
    assignments <- rep(0,numPoints);
    for(n in 1:numPoints) {
        max.cluster <- 0
        max.assignment <- -1
        # Prune any points that do not have an appreciable probability
        # of being in their respective cluster
        #max.assignment <- 1/numClusters
        #if(numClusters > 2) {
        #  max.assignment <- (4/3)/numClusters
        #}
        for(k in 1:numClusters) {
            if ( !is.na(probabilities[n,k]) & (probabilities[n,k] > max.assignment) ) {
                max.assignment <- probabilities[n,k]
                max.cluster <- k
            }
        }
        assignments[n] <- max.cluster
    }
    return(assignments)
}


##--------------------------------------------------------------------------
## Do clustering with bmm (beta mixture model)
##
clusterWithBmm <- function(vafs.merged, vafs, vars, refs, initialClusters=10, samples=1, plotIntermediateResults=0, verbose=TRUE, overlap.threshold=0.8) {
    library(bmm)

    initialClusters=initialClusters
    if(length(vafs[,1]) <= initialClusters){
      print(paste("ERROR: only",length(vafs[,1])," points 0 not enough points to cluster when using",initialClusters,"intialClusters. Provide more data or reduce your maximumClusters option"))
      return(list(NULL))
    }

    #replace any values of zero with a very small number to prevent errors
    # NB:  It seems OK to add .Machine$double.eps (without the factor of 100)
    # to 0 to avoid trouble, but subtracting from 1 is not sufficient.  Presumably,
    # there isn't enough accuracy to represent such a number and it is stored as 1.
    delta <- 100 * .Machine$double.eps
    vafs[which(vafs <= delta)] = delta

    #also replace any values of one to prevent errors
    vafs[which(vafs >= (1-delta))] = 1 - delta

    ## Initialize the hyperparameters of the Beta mixture model (bmm).
    hyperparams <- init.bmm.hyperparameters(vafs, initialClusters)

    ## Initialize the parameters of the bmm.
    params <- init.bmm.parameters(vafs, initialClusters, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0)

    ## Perform the clustering.
    ## Start with the provided number of clusters, but prune any with low probability
    bmm.results <- bmm.filter.clusters(vafs.merged, vafs, initialClusters, params$r, params$mu, params$alpha, params$nu, params$beta, params$c, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose, plotIntermediateResults=plotIntermediateResults, overlap.threshold=overlap.threshold)
    if(bmm.results$retVal != 0) {
        cat("WARNING: bmm failed to converge. No clusters assigned\n")
        return(NULL);
    }

    ##get the assignment of each point to a cluster
    probs = bmm.results$r
    numPoints = length(probs[,1])
    numClusters = length(probs[1,])
    clusters = hardClusterAssignments(numPoints,numClusters,probs);

    ## find confidence intervals around the means of the clusters
    intervals = bmm.narrowest.mean.interval.about.centers(bmm.results$mu, bmm.results$alpha, bmm.results$nu, bmm.results$beta, 0.68)
    means = intervals$centers
    if(verbose){
      print("Cluster Centers:");
      print(means);
    }
    lower = intervals$lb
    upper = intervals$ub

    if(verbose){
      print("Outliers:")
      print(bmm.results$outliers)
    }


    ## Generate (x,y) values of the posterior predictive density
    n <- 1000

    ## Don't evaluate at x=0 or x=1, which will blow up
    x <- seq(0, 1, 1/n)[2:n]

    ## create a num_dimensions x n matrix of y values
    n=n-1;
    y <- rep.int(0, n)
    y = t(matrix(rep(y,dim(vafs)[2]),ncol=dim(vafs)[2]))

    ##for each dimension
    yms = list()
    for (k in 1:numClusters) {
      yms[[k]] <- y
    }
    for(dim in 1:dim(vafs)[2]){
        ym <- matrix(data=0, nrow=numClusters, ncol=n)
        num.iterations <- 100
        for (k in 1:numClusters) {
            for (i in 1:n) {
                ## Evaluate posterior probability at x.
                ym[k,i] <- bmm.component.posterior.predictive.density(x[i], bmm.results$mu[dim,k], bmm.results$alpha[dim,k], bmm.results$nu[dim,k], bmm.results$beta[dim,k], bmm.results$E.pi[k], num.samples = num.iterations)
                yms[[k]][dim,i] <- ym[k,i]
                y[dim,i] <- y[dim,i] + ym[k,i]
            }
        }
        ##scale yvals between 0 and 1
        #for (k in 1:numClusters) {
        #  yms[[k]][dim,] <- yms[[k]][dim,]/max(yms[[k]][dim,])
        #}
        #y[dim,] = y[dim,]/max(y[dim,])
    }

    ##scale xvals between 1 and 100
    x = x*100

    #return a list of info
    return(list(
        cluster.assignments = clusters,
        cluster.probabilities= probs,
        cluster.means = means,
        cluster.upper = upper,
        cluster.lower = lower,
        fit.x = x,
        fit.y = y,
        individual.fits.y = yms))
}

##--------------------------------------------------------------------------
## Do clustering with binomial bmm (binomial bayesian mixture model)
##
clusterWithBinomialBmm <- function(vafs.merged, vafs, vars, refs, initialClusters=10, samples=1, plotIntermediateResults=0, verbose=TRUE, overlap.threshold=0.8) {
    library(bmm)

    initialClusters=initialClusters
    if(length(vafs[,1]) <= initialClusters){
      print(paste("ERROR: only",length(vafs[,1])," points 0 not enough points to cluster when using",initialClusters,"intialClusters. Provide more data or reduce your maximumClusters option"))
      return(list(NULL))
    }

    total.trials <- vars + refs

    ## Initialize the hyperparameters of the Binomial mixture model.
    hyperparams <- init.binomial.bmm.hyperparameters(vars, total.trials, initialClusters)

    ## Initialize the parameters of the bmm.
    params <- init.binomial.bmm.parameters(vars, total.trials, initialClusters, hyperparams$a0, hyperparams$b0, hyperparams$alpha0)

    ## Perform the clustering.
    ## Start with the provided number of clusters, but prune any with low probability
    bmm.results <- binomial.bmm.filter.clusters(vafs.merged, vafs, vars, total.trials, initialClusters, params$r, params$a, params$b, params$alpha, hyperparams$a0, hyperparams$b0, hyperparams$alpha0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose, plotIntermediateResults=plotIntermediateResults, overlap.threshold=overlap.threshold)
    if(bmm.results$retVal != 0) {
        cat("WARNING: bmm failed to converge. No clusters assigned\n")
        return(NULL);
    }

    ##get the assignment of each point to a cluster
    probs = bmm.results$r
    numPoints = length(probs[,1])
    numClusters = length(probs[1,])
    clusters = hardClusterAssignments(numPoints,numClusters,probs);

    ## find confidence intervals around the means of the clusters
    intervals = binomial.bmm.narrowest.mean.interval.about.centers(bmm.results$a, bmm.results$b, 0.68)
    means = intervals$centers
    if(verbose){
      print("Cluster Centers:");
      print(means);
    }
    lower = intervals$lb
    upper = intervals$ub

    if(verbose){
      print("Outliers:")
      print(bmm.results$outliers)
    }


    ## Generate (x,y) values of the posterior predictive density

    ## We will evaluate these densities at the median total.trials
    ## in order to uniquely define proportions from count data.
    eta <- round(median(total.trials))

    x <- seq(0, eta, 1)
    n <- length(x)

    y <- rep.int(0, n)
    y = t(matrix(rep(y,dim(vafs)[2]),ncol=dim(vafs)[2]))

    ##for each dimension
    yms = list()
    for (k in 1:numClusters) {
      yms[[k]] <- y
    }
    for(dim in 1:dim(vafs)[2]){
        ym <- matrix(data=0, nrow=numClusters, ncol=n)
        num.iterations <- 100
        for (k in 1:numClusters) {
            for (i in 1:n) {
                ## Evaluate posterior probability at x.
                ym[k,i] <- binomial.bmm.component.posterior.predictive.density(x[i], eta, bmm.results$a[k,dim], bmm.results$b[k,dim], bmm.results$E.pi[k])
                yms[[k]][dim,i] <- ym[k,i]
                y[dim,i] <- y[dim,i] + ym[k,i]
            }
        }
    }

    x = ( x / eta ) * 100

    return(list(
        cluster.assignments = clusters,
        cluster.probabilities= probs,
        cluster.means = means,
        cluster.upper = upper,
        cluster.lower = lower,
        fit.x = x,
        fit.y = y,
        individual.fits.y = yms))
}


##--------------------------------------------------------------------------
## Do clustering with gaussian bmm (gaussian bayesian mixture model)
##
clusterWithGaussianBmm <- function(vafs.merged, vafs, vars, refs, initialClusters=10, samples=1, plotIntermediateResults=0, verbose=TRUE, overlap.threshold=0.8) {
    library(bmm)

    initialClusters=initialClusters
    if(length(vafs[,1]) <= initialClusters){
      print(paste("ERROR: only",length(vafs[,1])," points 0 not enough points to cluster when using",initialClusters,"intialClusters. Provide more data or reduce your maximumClusters option"))
      return(list(NULL))
    }

    total.trials <- vars + refs

    ## Initialize the hyperparameters of the gaussian mixture model.
    hyperparams <- init.gaussian.bmm.hyperparameters(vafs, initialClusters)

    ## Initialize the parameters of the model.
    params <- init.gaussian.bmm.parameters(vafs, initialClusters, hyperparams$m0, hyperparams$alpha0, hyperparams$beta0, hyperparams$nu0, hyperparams$W0)

    ## Perform the clustering.
    ## Start with the provided number of clusters, but prune any with low probability
    bmm.results <- gaussian.bmm.filter.clusters(vafs.merged, vafs, vars, total.trials, initialClusters, params$r, params$m, params$alpha, params$beta, params$nu, params$W, hyperparams$m0, hyperparams$alpha0, hyperparams$beta0, hyperparams$nu0, hyperparams$W0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose, plotIntermediateResults=plotIntermediateResults, overlap.threshold=overlap.threshold)
    if(bmm.results$retVal != 0) {
        cat("WARNING: bmm failed to converge. No clusters assigned\n")
        return(NULL);
    }

    ##get the assignment of each point to a cluster
    probs = bmm.results$r
    numPoints = length(probs[,1])
    numClusters = length(probs[1,])
    clusters = hardClusterAssignments(numPoints,numClusters,probs);

    ## find confidence intervals around the means of the clusters
    intervals = gaussian.bmm.narrowest.mean.interval.about.centers(bmm.results$m, bmm.results$alpha, bmm.results$beta, bmm.results$nu, bmm.results$W, 0.68)
    means = intervals$centers
    if(verbose){
      print("Cluster Centers:");
      print(means);
    }
    lower = intervals$lb
    upper = intervals$ub

    if(verbose){
      print("Outliers:")
      print(bmm.results$outliers)
    }


    # Only generate the fit values in 1D--that's the only
    # place we use it.  Unlike the beta and binomial approaches,
    # the dimensions are not independent here, so we would either
    # have to generate the fits in the full dimensional space
    # or we would have to integrate out all but one dimension.
    # We should be able to do that integration analytically,
    # but I'm lazy, so won't attempt it now.
    x <- NULL
    y <- NULL
    yms <- NULL

    if(dim(vafs)[2]==1) {
      ## Generate (x,y) values of the posterior predictive density
      n <- 1000

      ## Don't evaluate at x=0 or x=1, which will blow up
      x <- seq(0, 1, 1/n)[2:n]

      ## create a num_dimensions x n matrix of y values
      n=n-1;
      y <- rep.int(0, n)
      y = t(matrix(rep(y,dim(vafs)[2]),ncol=dim(vafs)[2]))

      ##for each dimension
      yms = list()
      for (k in 1:numClusters) {
        yms[[k]] <- y
      }
      for(dim in 1:dim(vafs)[2]){
          ym <- matrix(data=0, nrow=numClusters, ncol=n)
          num.iterations <- 100
          for (k in 1:numClusters) {
              for (i in 1:n) {
                  ## Evaluate posterior probability at x.
                  ym[k,i] <- gaussian.bmm.component.posterior.predictive.density(x[i], k, bmm.results$m, bmm.results$alpha, bmm.results$beta, bmm.results$nu, bmm.results$W, bmm.results$L)
                  yms[[k]][dim,i] <- ym[k,i]
                  y[dim,i] <- y[dim,i] + ym[k,i]
              }
          }
          ##scale yvals between 0 and 1
          #for (k in 1:numClusters) {
          #  yms[[k]][dim,] <- yms[[k]][dim,]/max(yms[[k]][dim,])
          #}
          #y[dim,] = y[dim,]/max(y[dim,])
      }

      ##scale xvals between 1 and 100
      x = x*100
    }

    #return a list of info
    return(list(
        cluster.assignments = clusters,
        cluster.probabilities= probs,
        cluster.means = means,
        cluster.upper = upper,
        cluster.lower = lower,
        fit.x = x,
        fit.y = y,
        individual.fits.y = yms))
}

## ##--------------------------------------------------------------------------
## ## The beta distribution clustering + filtering method
## ##
bmm.filter.clusters <- function(vafs.merged, X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0,
                                convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0,
                                plotIntermediateResults = 0, overlap.threshold=0.8){

    total.iterations <- 0
    N <- dim(X)[1]
    num.dimensions <- dim(X)[2]

    x.colnames <- colnames(X)

    outliers <- matrix(data=0, nrow=0, ncol=dim(X)[2])
    colnames(outliers) <- colnames(X)

    E.pi.prev <- rep(0, N.c)

    width <- as.double(erf(1.5/sqrt(2)))
    # width <- as.double(erf(1/sqrt(2)))
    if(plotIntermediateResults > 0) {

      probs <- r
      numPoints = length(probs[,1])
      numClusters = length(probs[1,])
      clusters = hardClusterAssignments(numPoints,numClusters,probs);
      vafs.with.assignments = cbind(vafs.merged,cluster=clusters)
      outputPrefix <- paste("tumor", total.iterations, sep="")
      sampleNames <- c("Tumor", "Relapse")
      positionsToHighlight <- NULL
      highlightsHaveNames <- FALSE
      overlayClusters <- TRUE

      ellipse.width <- as.double(erf(1/sqrt(2)))

      # Calculate standard error of the means
      SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, ellipse.width)
      SEM.centers <- 100 * t(SEM.res$centers)
      SEMs.lb <- 100 * t(SEM.res$lb)
      SEMs.ub <- 100 * t(SEM.res$ub)

      # Calculate standard errors
      std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, ellipse.width)
      std.dev.centers <- 100 * t(std.dev.res$centers)
      std.dev.lb <- 100 * t(std.dev.res$lb)
      std.dev.ub <- 100 * t(std.dev.res$ub)

      ellipse.metadata <- list(SEMs.lb = SEMs.lb, SEMs.ub = SEMs.ub, std.dev.lb = std.dev.lb, std.dev.ub = std.dev.ub)

      #sc.plot2d(vafs.with.assignments, outputPrefix, sampleNames, dim(X)[2], ellipse.metadata=ellipse.metadata)
      sc.plot2d(vafs.with.assignments, outputPrefix, ellipse.metadata=ellipse.metadata, positionsToHighlight=positionsToHighlight, highlightsHaveNames=highlightsHaveNames)
    }

    # Outer while loop: following convergence of inner loop, apply
    # overlapping cluster condition to drop any overlapping clusters.
    while(TRUE) {

        #print(dim(X))
        #print(N.c)
        #print(N)
        #print(length(c0))
        #print(dim(r))
        #print(length(colSums(r)))

        if(plotIntermediateResults > 0) {
          max.iterations <- plotIntermediateResults
        }

        bmm.res <- bmm.fixed.num.components(X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0, convergence.threshold, max.iterations, verbose)
        if((bmm.res$retVal != 0) & (plotIntermediateResults == 0)) {
            cat("Failed to converge!\n")
            q(status=-1)
        }
        mu <- bmm.res$mu
        alpha <- bmm.res$alpha
        nu <- bmm.res$nu
        beta <- bmm.res$beta
        c <- bmm.res$c
        E.pi <- bmm.res$E.pi
        ubar <- bmm.res$ubar
        vbar <- bmm.res$vbar
        r <- bmm.res$r

        total.iterations <- total.iterations + bmm.res$num.iterations

        ln.rho <- bmm.res$ln.rho
        E.lnu <- bmm.res$E.lnu
        E.lnv <- bmm.res$E.lnv
        E.lnpi <- bmm.res$E.lnpi
        E.quadratic.u <- bmm.res$E.quadratic.u
        E.quadratic.v <- bmm.res$E.quadratic.v

        if(plotIntermediateResults > 0) {
          probs <- r
          numPoints = length(probs[,1])
          numClusters = length(probs[1,])
          clusters = hardClusterAssignments(numPoints,numClusters,probs);
          vafs.with.assignments = cbind(vafs.merged,cluster=clusters)
          #outputPrefix <- paste("tumor", total.iterations, sep="")
          sampleNames <- c("Tumor", "Relapse")
          positionsToHighlight <- NULL
          highlightsHaveNames <- FALSE
          overlayClusters <- TRUE

          ellipse.width <- as.double(erf(1/sqrt(2)))

          # Calculate standard error of the means
          SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, ellipse.width)
          SEM.centers <- 100 * t(SEM.res$centers)
          SEMs.lb <- 100 * t(SEM.res$lb)
          SEMs.ub <- 100 * t(SEM.res$ub)

          # Calculate standard errors
          std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, ellipse.width)
          std.dev.centers <- 100 * t(std.dev.res$centers)
          std.dev.lb <- 100 * t(std.dev.res$lb)
          std.dev.ub <- 100 * t(std.dev.res$ub)

          ellipse.metadata <- list(SEMs.lb = SEMs.lb, SEMs.ub = SEMs.ub, std.dev.lb = std.dev.lb, std.dev.ub = std.dev.ub)
          ##todo - fixme
          ##sc.plot2d(vafs.with.assignments, "clusters.int.pdf", sampleNames, dim(X)[2], positionsToHighlight, highlightsHaveNames, overlayClusters, ellipse.metadata=ellipse.metadata)

        }

        if((bmm.res$retVal != 0) & (plotIntermediateResults > 0)) {
          next
        }


        do.inner.iteration <- FALSE

        # Remove any small clusters

        apply.min.items.condition <- TRUE
        apply.uncertainty.self.overlap.condition <- FALSE
        apply.large.SEM.condition <- FALSE
        apply.overlapping.SEM.condition <- FALSE
        apply.overlapping.std.dev.condition <- TRUE

        # Changes on Mar 25, 2013
        apply.overlapping.std.dev.condition <- FALSE
        apply.uncertainty.self.overlap.condition <- TRUE

        if((apply.min.items.condition == TRUE) & (N.c > 1)) {
            # threshold.pts <- 10
            threshold.pts <- 3

            clusters <- hardClusterAssignments(N,N.c,r)

            num.items.per.cluster <- rep(0, N.c)
            for(n in 1:N) {
                num.items.per.cluster[clusters[n]] <- num.items.per.cluster[clusters[n]] + 1
            }

            non.zero.indices <- num.items.per.cluster >= threshold.pts

            if ( any(non.zero.indices==FALSE) ) {

                do.inner.iteration <- TRUE

                numeric.indices <- (1:N.c)
                if(verbose){
                  cat("Dropping clusters with pts: \n")
                  print(num.items.per.cluster[!non.zero.indices])
                }

                remove.data <- TRUE
                if (remove.data == TRUE) {
                    numeric.indices <- (1:N.c)[non.zero.indices]
                    new.outliers <- matrix(X[!(clusters %in% numeric.indices),], ncol=num.dimensions)
                    outliers <- rbind(outliers, new.outliers)
                }

                # To remove data from the data set, set its entries to NA
                # X <- matrix(X[clusters %in% numeric.indices,], ncol=num.dimensions)
                # N <- dim(X)[1]

                X[!(clusters %in% numeric.indices),] <- NA
                # N <- dim(X)[1] - dim(outliers)[1]

                colnames(X) <- x.colnames

                E.pi <- E.pi[non.zero.indices]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[non.zero.indices]
                c <- c[non.zero.indices]
                c0 <- c0[non.zero.indices]

                # debugging for the issue of mismatched matrix size #
                #print("line 248");
                #print(dim(r));
                #print(dim(r[,non.zero.indices]));
                #print(N);
                #print(N.c);
                #print(summary(factor(non.zero.indices)));

                # Don't resize r and ln.rho matrices to accomodate removed
                # outliers, instead we will have set their rows to NA above.
                # r <- matrix(r[clusters %in% numeric.indices,non.zero.indices], nrow=N, ncol=N.c)
                # ln.rho <- matrix(ln.rho[clusters %in% numeric.indices,non.zero.indices], nrow=N, ncol=N.c)
                # But do drop any columns corresponding to dropped clusters
                r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
                    if(any(is.na(ln.rho[n,]))) {
                      r[n,] <- rep(NA, N.c)
                      next
                    }

                    row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
                    for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
                }

                mu <- matrix(mu[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                nu <- matrix(nu[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                mu0 <- matrix(mu0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                nu0 <- matrix(nu0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                alpha <- matrix(alpha[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                beta <- matrix(beta[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                alpha0 <- matrix(alpha0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                beta0 <- matrix(beta0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                ubar <- matrix(ubar[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                vbar <- matrix(vbar[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                E.pi.prev <- E.pi

                E.lnpi <- E.lnpi[non.zero.indices]
                E.lnu <- E.lnu[non.zero.indices]
                E.lnv <- E.lnv[non.zero.indices]
                E.quadratic.u <- matrix(E.quadratic.u[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                E.quadratic.v <- matrix(E.quadratic.v[,non.zero.indices], nrow=num.dimensions, ncol=N.c)

            }
        } # End apply.min.items.condition

        if((apply.uncertainty.self.overlap.condition == TRUE) & (N.c > 1)) {

            overlaps <- rep(0, N.c)
            ones <- rep(1, N)
            indices.to.keep.boolean <- rep(TRUE, N.c)
            # Just drop min overlap
            for(k in 1:N.c) {
                overlaps[k] <- sum(r[,k] * r[,k], na.rm=TRUE) / sum(ones * r[,k], na.rm=TRUE)
            }

            if(min(overlaps, na.rm=TRUE) < overlap.threshold) {
                for(k in 1:N.c) {
                    # Small clusters will have undefined overlaps, just skip.
                    # We'll remove them later.
                    if(is.nan(overlaps[k])) { next }
                    if((overlaps[k] < overlap.threshold) & (overlaps[k] == min(overlaps, na.rm=TRUE))) {
                        indices.to.keep.boolean[k] <- FALSE
                        break
                    }
                }
            }

            if(verbose){
              for(k in 1:N.c) {
                cat(sprintf("Cluster %d pi = %.3f self-overlap = %.3f\n", k, E.pi[k], overlaps[k]))
              }
            }
            indices.to.keep <- (1:N.c)[indices.to.keep.boolean]

            if(length(indices.to.keep) != N.c) {

                do.inner.iteration <- TRUE

                means <- matrix(ubar/(ubar+vbar), nrow=num.dimensions, ncol=N.c)
                if(length(indices.to.keep) < N.c) {
                    indices.to.drop <- (1:N.c)[!indices.to.keep.boolean]
                    for(i in 1:length(indices.to.drop)) {
                        index <- indices.to.drop[i]
                        if(verbose){
                          cat("Self-overlap condition: Dropping cluster with center: ")
                          for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                          }
                          cat("\n")
                        }
                    }
                }

                # NB: only removing clusters here, not data points
                E.pi <- E.pi[indices.to.keep]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[indices.to.keep]
                c <- c[indices.to.keep]
                c0 <- c0[indices.to.keep]

                r <- matrix(r[,indices.to.keep], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,indices.to.keep], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
                    if(any(is.na(ln.rho[n,]))) {
                      r[n,] <- rep(NA, N.c)
                      next
                    }

                    row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
                    for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
                }

                mu <- matrix(mu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu <- matrix(nu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                mu0 <- matrix(mu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu0 <- matrix(nu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha <- matrix(alpha[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta <- matrix(beta[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha0 <- matrix(alpha0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta0 <- matrix(beta0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                ubar <- matrix(ubar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                vbar <- matrix(vbar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.pi.prev <- E.pi

                E.lnpi <- E.lnpi[indices.to.keep]
                E.lnu <- E.lnu[indices.to.keep]
                E.lnv <- E.lnv[indices.to.keep]
                E.quadratic.u <- matrix(E.quadratic.u[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.quadratic.v <- matrix(E.quadratic.v[,indices.to.keep], nrow=num.dimensions, ncol=N.c)


            }
        } # End apply.uncertainty.self.overlap.condition

        if((apply.large.SEM.condition == TRUE) & (N.c > 1)) {

            # Calculate standard error of the means
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)

            # Calculate standard errors
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)

            indices.to.keep.boolean <- rep(TRUE, N.c)

            if(verbose){
              show.output <- TRUE
            }
            for(k in 1:N.c) {
                if (show.output) {
                    cat(sprintf("%dD: Cluster %d pi = %.3f: ", num.dimensions, k, E.pi[k]))
                }
                greater.than.30 <- TRUE
                greater.than.02 <- TRUE
                SEM.width.threshold <- 0.02
                #SEM.width.threshold <- 0.001
                for(m in 1:num.dimensions) {
                    center <- SEM.centers[k,m]
                    lower <- SEMs.lb[k,m]
                    upper <- SEMs.ub[k,m]
                    std.dev.width <- (std.dev.ub[k,m] - std.dev.lb[k,m])
                    SEM.width <- (SEMs.ub[k,m] - SEMs.lb[k,m])
                    #if((SEM.width/std.dev.width)>.3){ greater.than.30 <- TRUE }
                    #if(SEM.width>.02){ greater.than.02 <- TRUE }
                    if((SEM.width/std.dev.width)<=.3){ greater.than.30 <- FALSE }
                    if(SEM.width<=SEM.width.threshold){ greater.than.02 <- FALSE }
                    if (show.output) {
                        cat(sprintf("%.3f (%.3f, %.3f) [(%.3f) %.3f, %.3f] {%.3f} ", center, lower, upper, width, std.dev.lb[k,m], std.dev.ub[k,m], SEM.width/std.dev.width))
                        if(greater.than.30) { cat("* ") }
                        if(greater.than.02) { cat("**") }
                    }
                }
                if (show.output) {
                    cat("\n")
                }
                if(greater.than.30 & greater.than.02) { indices.to.keep.boolean[k] <- FALSE }
            }

            non.zero.indices <- indices.to.keep.boolean

            if ( any(non.zero.indices==FALSE) ) {

                do.inner.iteration <- TRUE
                numeric.indices <- (1:N.c)

                if(verbose){
                  cat("Removing clusters because of SEM condition: \n")
                  print(numeric.indices[non.zero.indices])
                }

                clusters <- hardClusterAssignments(N,N.c,r)

                remove.data <- TRUE
                if (remove.data == TRUE) {
                    numeric.indices <- (1:N.c)[non.zero.indices]
                    new.outliers <- matrix(X[!(clusters %in% numeric.indices),], ncol=num.dimensions)
                    outliers <- rbind(outliers, new.outliers)
                }

                # To remove data from the data set, set its entries to NA
                # X <- matrix(X[clusters %in% numeric.indices,], ncol=num.dimensions)
                # N <- dim(X)[1]
                X[!(clusters %in% numeric.indices),] <- NA
                # N <- dim(X)[1] - dim(outliers)[1]

                colnames(X) <- x.colnames

                E.pi <- E.pi[non.zero.indices]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[non.zero.indices]
                c <- c[non.zero.indices]
                c0 <- c0[non.zero.indices]

                # Don't resize r and ln.rho matrices to accomodate removed
                # outliers, instead we will have set their rows to NA above.
                # r <- matrix(r[clusters %in% numeric.indices,non.zero.indices], nrow=N, ncol=N.c)
                # ln.rho <- matrix(ln.rho[clusters %in% numeric.indices,non.zero.indices], nrow=N, ncol=N.c)
                # But do drop any columns corresponding to dropped clusters
                r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)

                # Need to renormalize r--do it gently.
                for(n in 1:N) {
                    if(any(is.na(ln.rho[n,]))) {
                      r[n,] <- rep(NA, N.c)
                      next
                    }

                    row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
                    for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
                }

                mu <- matrix(mu[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                nu <- matrix(nu[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                mu0 <- matrix(mu0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                nu0 <- matrix(nu0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                alpha <- matrix(alpha[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                beta <- matrix(beta[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                alpha0 <- matrix(alpha0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                beta0 <- matrix(beta0[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                ubar <- matrix(ubar[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                vbar <- matrix(vbar[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                E.pi.prev <- E.pi

                E.lnpi <- E.lnpi[non.zero.indices]
                E.lnu <- E.lnu[non.zero.indices]
                E.lnv <- E.lnv[non.zero.indices]
                E.quadratic.u <- matrix(E.quadratic.u[,non.zero.indices], nrow=num.dimensions, ncol=N.c)
                E.quadratic.v <- matrix(E.quadratic.v[,non.zero.indices], nrow=num.dimensions, ncol=N.c)


            }

        } # End apply.large.SEM.condition

        if((apply.overlapping.SEM.condition == TRUE) & (N.c > 1)) {

            # Calculate standard error of the means
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)

            # Calculate standard errors
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)

            indices.to.keep.boolean <- rep(TRUE, N.c)
            # Determine if component i's center is contained within
            # component i2's std.dev
            pi.threshold = 10^-2

            for(i in 1:N.c){
                i.subsumed.by.another.cluster <- FALSE
                for(i2 in 1:N.c){
                    if(i == i2) { next }
                    if(indices.to.keep.boolean[i2] == FALSE) { next }
                    if(E.pi[i2] < pi.threshold) { next }
                    i.subsumed.by.another.cluster <- TRUE
                    for(l in 1:num.dimensions){
                        i.center <- std.dev.centers[i,l]
                        if((i.center < std.dev.lb[i2,l]) | (i.center > std.dev.ub[i2,l])) {
                            i.subsumed.by.another.cluster <- FALSE
                            break
                        }
                    }
                    if(i.subsumed.by.another.cluster==TRUE) {
                      if(verbose){
                        cat(sprintf("2. Dropping cluster with center: "))
                        for(l in 1:num.dimensions){
                          cat(sprintf("%.3f ", std.dev.centers[i,l]))
                        }
                        cat(sprintf("because it overlaps with: "))
                        for(l in 1:num.dimensions){
                          cat(sprintf("(%.3f, %.3f) ", std.dev.lb[i2,l], std.dev.ub[i2,l]))
                        }
                        cat("\n")
                        break
                      }
                    }
                }
                if(i.subsumed.by.another.cluster==TRUE) {
                    indices.to.keep.boolean[i] <- FALSE
                }
            }

            indices.to.keep <- (1:N.c)[indices.to.keep.boolean]

            if(length(indices.to.keep) != N.c) {

                do.inner.iteration <- TRUE

                means <- matrix(ubar/(ubar+vbar), nrow=num.dimensions, ncol=N.c)
                if(length(indices.to.keep) < N.c) {
                    indices.to.drop <- (1:N.c)[!indices.to.keep.boolean]
                    for(i in 1:length(indices.to.drop)) {
                        index <- indices.to.drop[i]
                        if(verbose){
                          cat("3. Dropping cluster with center: ")
                          for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                          }
                          cat("\n")
                        }
                    }
                }

                # NB: only removing clusters here, not data points
                E.pi <- E.pi[indices.to.keep]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[indices.to.keep]
                c <- c[indices.to.keep]
                c0 <- c0[indices.to.keep]

                r <- matrix(r[,indices.to.keep], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,indices.to.keep], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
                    if(any(is.na(ln.rho[n,]))) {
                      r[n,] <- rep(NA, N.c)
                      next
                    }

                    row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
                    for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
                }

                mu <- matrix(mu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu <- matrix(nu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                mu0 <- matrix(mu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu0 <- matrix(nu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha <- matrix(alpha[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta <- matrix(beta[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha0 <- matrix(alpha0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta0 <- matrix(beta0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                ubar <- matrix(ubar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                vbar <- matrix(vbar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.pi.prev <- E.pi

                E.lnpi <- E.lnpi[indices.to.keep]
                E.lnu <- E.lnu[indices.to.keep]
                E.lnv <- E.lnv[indices.to.keep]
                E.quadratic.u <- matrix(E.quadratic.u[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.quadratic.v <- matrix(E.quadratic.v[,indices.to.keep], nrow=num.dimensions, ncol=N.c)

            }

        } # End apply.overlapping.SEM.condition

        if((apply.overlapping.std.dev.condition == TRUE) & (N.c > 1)) {

            # Calculate standard error of the means
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)

            # Calculate standard errors
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)

            indices.to.keep.boolean <- rep(TRUE, N.c)
            # Determine how much component i's std dev overlaps component i2's std dev
            overlaps <- matrix(data = 1, nrow=N.c, ncol=N.c)
            std.dev.overlap.threshold <- 0
            for(i in 1:N.c){
                #cat("Cluster i = ", i, " overlaps: \n")
                for(i2 in 1:N.c){
                  if(verbose){
                    cat("   ")
                  }
                    for(l in 1:num.dimensions){
                        fraction <- 0
                        if((std.dev.lb[i,l] < std.dev.ub[i2,l]) & (std.dev.ub[i,l] > std.dev.lb[i2,l])) {
                            overlap <- min(std.dev.ub[i,l], std.dev.ub[i2,l]) - max(std.dev.lb[i,l], std.dev.lb[i2,l])
                            fraction <- overlap / ( std.dev.ub[i,l] - std.dev.lb[i,l] )
                        }
                        if((fraction > std.dev.overlap.threshold) & (overlaps[i,i2] != 0)) {
                            overlaps[i,i2] <- 1
                        } else {
                            overlaps[i,i2] <- 0
                        }
                        #cat(sprintf("f = %.3f ", fraction))
                    }
                    #cat("\n")
                }
            }

            # Remove one of two overlapping clusters.  If both overlap
            # (NB: above is not symmetric), then only remove smaller of two.
            for(i in 2:N.c){
                if(indices.to.keep.boolean[i] == FALSE) { next }
                for(i2 in 1:(i-1)){
                    if(indices.to.keep.boolean[i2] == FALSE) { next }
                    if(overlaps[i,i2] == 1) {
                        if((overlaps[i2,i] == 1) & (E.pi[i2] < E.pi[i])) {
                          if(verbose){cat("Removing ", i2, " because of overlap with i = ", i, "\n")}
                            indices.to.keep.boolean[i2] <- FALSE
                        } else {
                            if(verbose){cat("Removing ", i, " because of overlap with i2 = ", i2, "\n")}
                            indices.to.keep.boolean[i] <- FALSE
                        }
                    }
                }
            }

            indices.to.keep <- (1:N.c)[indices.to.keep.boolean]

            if(length(indices.to.keep) != N.c) {

                do.inner.iteration <- TRUE

                means <- matrix(ubar/(ubar+vbar), nrow=num.dimensions, ncol=N.c)
                if(length(indices.to.keep) < N.c) {
                    indices.to.drop <- (1:N.c)[!indices.to.keep.boolean]
                    for(i in 1:length(indices.to.drop)) {
                        index <- indices.to.drop[i]
                        if(verbose){
                          cat("4. Dropping cluster with center: ")
                          for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                          }
                          cat("\n")
                        }
                    }
                }

                # NB: only removing clusters here, not data points
                E.pi <- E.pi[indices.to.keep]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[indices.to.keep]
                c <- c[indices.to.keep]
                c0 <- c0[indices.to.keep]

                r <- matrix(r[,indices.to.keep], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,indices.to.keep], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
                    if(any(is.na(ln.rho[n,]))) {
                      r[n,] <- rep(NA, N.c)
                      next
                    }

                    row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
                    for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
                }

                mu <- matrix(mu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu <- matrix(nu[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                mu0 <- matrix(mu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                nu0 <- matrix(nu0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha <- matrix(alpha[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta <- matrix(beta[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                alpha0 <- matrix(alpha0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                beta0 <- matrix(beta0[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                ubar <- matrix(ubar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                vbar <- matrix(vbar[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.pi.prev <- E.pi

                E.lnpi <- E.lnpi[indices.to.keep]
                E.lnu <- E.lnu[indices.to.keep]
                E.lnv <- E.lnv[indices.to.keep]
                E.quadratic.u <- matrix(E.quadratic.u[,indices.to.keep], nrow=num.dimensions, ncol=N.c)
                E.quadratic.v <- matrix(E.quadratic.v[,indices.to.keep], nrow=num.dimensions, ncol=N.c)

            }

        } # End apply.overlapping.std.dev.condition


        if(do.inner.iteration == FALSE) { break }

    } # End outer while(TRUE)

    # Calculate standard error of the means
    SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, width)
    SEM.centers <- t(SEM.res$centers)
    SEMs.lb <- t(SEM.res$lb)
    SEMs.ub <- t(SEM.res$ub)

    # Calculate standard errors
    std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, width)
    std.dev.centers <- t(std.dev.res$centers)
    std.dev.lb <- t(std.dev.res$lb)
    std.dev.ub <- t(std.dev.res$ub)

    for(k in 1:N.c) {
      if(verbose){cat(sprintf("Cluster %d pi = %.3f center =", k, E.pi[k]))}

        for(d in 1:num.dimensions) {
            center <- SEM.centers[k,d]
            if(verbose){cat(sprintf(" %.3f", center))}
        }
        if(verbose){cat(" SEM =")}
        for(d in 1:num.dimensions) {
            lb <- SEMs.lb[k,d]
            ub <- SEMs.ub[k,d]
            if(verbose){cat(sprintf(" (%.3f, %.3f)", lb, ub))}
        }
        if(verbose){cat(" sd =")}
        for(d in 1:num.dimensions) {
            lb <- std.dev.lb[k,d]
            ub <- std.dev.ub[k,d]
            if(verbose){cat(sprintf(" (%.3f, %.3f)", lb, ub))}
        }
        if(verbose){cat("\n")}
    }

    if(verbose){cat(sprintf('total iterations = %d\n', total.iterations))}

    retList <- list("retVal" = 0, "mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnu" = E.lnu, "E.lnv" = E.lnv, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "E.quadratic.u" = E.quadratic.u, "E.quadratic.v" = E.quadratic.v, "ubar" = ubar, "vbar" = vbar, "outliers" = outliers)

    return(retList)

} # End bmm.filter.clusters



## ##--------------------------------------------------------------------------
## ## The binomial distribution clustering + filtering method
## ##
binomial.bmm.filter.clusters <- function(vafs.merged, vafs, successes, total.trials, N.c, r, a, b, alpha, a0, b0, alpha0,
                                convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0,
                                plotIntermediateResults = 0, overlap.threshold=0.8){


  total.iterations <- 0
  num.dimensions <- dim(successes)[2]

  N <- dim(successes)[1]

  successes.colnames <- colnames(successes)
  total.colnames <- colnames(total.trials)

  outliers <- matrix(data=0, nrow=0, ncol=num.dimensions)
  colnames(outliers) <- colnames(vafs)

  E.pi.prev <- rep(0, N.c)

  width <- as.double(erf(1.5/sqrt(2)))

  while(TRUE) {

    if(verbose){
      print(r)
    }

    bmm.res <- binomial.bmm.fixed.num.components(successes, total.trials, N.c, r, a, b, alpha, a0, b0, alpha0, convergence.threshold, max.iterations = 10000, verbose = verbose)
    if(bmm.res$retVal != 0) {
      cat("Failed to converge!\n")
      q(status=-1)
    }

    a <- bmm.res$a
    b <- bmm.res$b
    alpha <- bmm.res$alpha

    E.pi <- bmm.res$E.pi
    r <- bmm.res$r

    total.iterations <- total.iterations + bmm.res$num.iterations

    ln.rho <- bmm.res$ln.rho
    E.lnpi <- bmm.res$E.lnpi

    do.inner.iteration <- FALSE

    apply.min.items.condition <- TRUE
    apply.uncertainty.self.overlap.condition <- FALSE
    apply.large.SEM.condition <- FALSE
    apply.overlapping.SEM.condition <- FALSE
    apply.overlapping.std.dev.condition <- TRUE

    # Changes on Mar 25, 2013
    apply.overlapping.std.dev.condition <- FALSE
    apply.uncertainty.self.overlap.condition <- TRUE

    # remove.data = TRUE iff we should remove points assigned to clusters
    # that we remove
    remove.data <- FALSE

    clusters <- hardClusterAssignments(N,N.c,r)

    indices.to.keep <- rep(TRUE, N.c)
    if((apply.min.items.condition == TRUE) & (N.c > 1)) {
      threshold.pts <- 10

      num.items.per.cluster <- rep(0, N.c)
      for(n in 1:N) {
        num.items.per.cluster[clusters[n]] <- num.items.per.cluster[clusters[n]] + 1
      }

      indices.to.keep <- num.items.per.cluster >= threshold.pts

      # indices.to.keep <- E.pi > pi.threshold
    } # End apply.min.items.condition

    if(all(indices.to.keep==TRUE) & (apply.uncertainty.self.overlap.condition == TRUE) & (N.c > 1)) {

      overlaps <- rep(0, N.c)
      ones <- rep(1, N)
      # Just drop min overlap
      for(k in 1:N.c) {
        overlaps[k] <- sum(r[,k] * r[,k], na.rm=TRUE) / sum(ones * r[,k], na.rm=TRUE)
      }

      if(min(overlaps, na.rm=TRUE) < overlap.threshold) {
        for(k in 1:N.c) {
          # Small clusters will have undefined overlaps, just skip.
          # We'll remove them later.
          if(is.nan(overlaps[k])) { next }
          if((overlaps[k] < overlap.threshold) & (overlaps[k] == min(overlaps, na.rm=TRUE))) {
            indices.to.keep[k] <- FALSE
            break
          }
        }
      }

      if(verbose){
        for(k in 1:N.c) {
          cat(sprintf("Cluster %d pi = %.3f self-overlap = %.3f\n", k, E.pi[k], overlaps[k]))
        }

        means <- matrix(a/(a+b), nrow=num.dimensions, ncol=N.c)
        indices.to.drop <- (1:N.c)[!indices.to.keep]
        for(i in 1:length(indices.to.drop)) {
          index <- indices.to.drop[i]
          if(verbose){
            cat("Self-overlap condition: Dropping cluster with center: ")
            for(l in 1:num.dimensions) {
              cat(sprintf("%.3f ", means[l, index]))
            }
            cat("\n")
          }
        }
      }

    } # End apply.uncertainty.self.overlap.condition

    # Calculate standard error of the means
    SEM.res <- binomial.bmm.narrowest.mean.interval.about.centers(a, b, width)
    SEM.centers <- t(SEM.res$centers)
    SEMs.lb <- t(SEM.res$lb)
    SEMs.ub <- t(SEM.res$ub)

    if(all(indices.to.keep==TRUE) & (apply.large.SEM.condition == TRUE) & (N.c > 1)) {

      for(k in 1:N.c) {
        if (verbose) {
          cat(sprintf("%dD: Cluster %d pi = %.3f: ", num.dimensions, k, E.pi[k]))
        }
        greater.than.30 <- TRUE
        greater.than.02 <- TRUE
        SEM.width.threshold <- 0.02
        #SEM.width.threshold <- 0.001
        for(m in 1:num.dimensions) {
          center <- SEM.centers[k,m]
          lower <- SEMs.lb[k,m]
          upper <- SEMs.ub[k,m]
          std.dev.width <- (std.dev.ub[k,m] - std.dev.lb[k,m])
          SEM.width <- (SEMs.ub[k,m] - SEMs.lb[k,m])
          #if((SEM.width/std.dev.width)>.3){ greater.than.30 <- TRUE }
          #if(SEM.width>.02){ greater.than.02 <- TRUE }
          if((SEM.width/std.dev.width)<=.3){ greater.than.30 <- FALSE }
          if(SEM.width<=SEM.width.threshold){ greater.than.02 <- FALSE }
          if (verbose) {
            cat(sprintf("%.3f (%.3f, %.3f) [(%.3f) %.3f, %.3f] {%.3f} ", center, lower, upper, width, std.dev.lb[k,m], std.dev.ub[k,m], SEM.width/std.dev.width))
            if(greater.than.30) { cat("* ") }
            if(greater.than.02) { cat("**") }
          }
        }
        if (verbose) {
          cat("\n")
        }
        if(greater.than.30 & greater.than.02) { indices.to.keep[k] <- FALSE }
      }
      if ( any(indices.to.keep==FALSE) ) {
        remove.data <- TRUE
      }
    } # End apply.large.SEM.condition

    if(all(indices.to.keep==TRUE) & (apply.overlapping.SEM.condition == TRUE) & (N.c > 1)) {

      # Determine if component i's center is contained within
      # component i2's std.dev
      pi.threshold = 10^-2

      for(i in 1:N.c){
        i.subsumed.by.another.cluster <- FALSE
        for(i2 in 1:N.c){
          if(i == i2) { next }
          if(indices.to.keep[i2] == FALSE) { next }
          if(E.pi[i2] < pi.threshold) { next }
          i.subsumed.by.another.cluster <- TRUE
          for(l in 1:num.dimensions){
            i.center <- std.dev.centers[i,l]
            if((i.center < std.dev.lb[i2,l]) | (i.center > std.dev.ub[i2,l])) {
              i.subsumed.by.another.cluster <- FALSE
              break
            }
          }
          if(i.subsumed.by.another.cluster==TRUE) {
            if(verbose){
              cat(sprintf("2. Dropping cluster with center: "))
              for(l in 1:num.dimensions){
                cat(sprintf("%.3f ", std.dev.centers[i,l]))
              }
              cat(sprintf("because it overlaps with: "))
              for(l in 1:num.dimensions){
                cat(sprintf("(%.3f, %.3f) ", std.dev.lb[i2,l], std.dev.ub[i2,l]))
              }
              cat("\n")
              break
            }
          }
        }
        if(i.subsumed.by.another.cluster==TRUE) {
          indices.to.keep[i] <- FALSE
        }
      }

    } # End apply.overlapping.SEM.condition

    if(all(indices.to.keep==TRUE) & (apply.overlapping.std.dev.condition == TRUE) & (N.c > 1)) {

      # Determine how much component i's std dev overlaps component i2's std dev
      overlaps <- matrix(data = 1, nrow=N.c, ncol=N.c)
      std.dev.overlap.threshold <- 0
      for(i in 1:N.c){
        for(i2 in 1:N.c){
          if(verbose){
            cat("   ")
          }
          for(l in 1:num.dimensions){
            fraction <- 0
            if((std.dev.lb[i,l] < std.dev.ub[i2,l]) & (std.dev.ub[i,l] > std.dev.lb[i2,l])) {
              overlap <- min(std.dev.ub[i,l], std.dev.ub[i2,l]) - max(std.dev.lb[i,l], std.dev.lb[i2,l])
              fraction <- overlap / ( std.dev.ub[i,l] - std.dev.lb[i,l] )
            }
            if((fraction > std.dev.overlap.threshold) & (overlaps[i,i2] != 0)) {
              overlaps[i,i2] <- 1
            } else {
              overlaps[i,i2] <- 0
            }
          }
        }
      }

      # Remove one of two overlapping clusters.  If both overlap
      # (NB: above is not symmetric), then only remove smaller of two.
      for(i in 2:N.c){
        if(indices.to.keep[i] == FALSE) { next }
        for(i2 in 1:(i-1)){
          if(indices.to.keep[i2] == FALSE) { next }
          if(overlaps[i,i2] == 1) {
            if((overlaps[i2,i] == 1) & (E.pi[i2] < E.pi[i])) {
              if(verbose){cat("Removing ", i2, " because of overlap with i = ", i, "\n")}
              indices.to.keep[i2] <- FALSE
            } else {
              if(verbose){cat("Removing ", i, " because of overlap with i2 = ", i2, "\n")}
              indices.to.keep[i] <- FALSE
            }
          }
        }
      }
    } # End apply.overlapping.std.dev.condition



    if ( any(indices.to.keep==FALSE) ) {

      do.inner.iteration <- TRUE

      numeric.indices <- (1:N.c)

      E.pi <- E.pi[indices.to.keep]
      E.lnpi <- E.lnpi[indices.to.keep]

      N.c <- length(E.pi)

      if(remove.data == TRUE) {
        cluster.indices.to.keep <- (1:N.c)[indices.to.keep]
        new.outliers <- matrix(vafs[!(clusters %in% cluster.indices.to.keep),], ncol=num.dimensions)
        outliers <- rbind(outliers, new.outliers)
        # To remove data from the data set, set its entries to NA
        vafs[!(clusters %in% cluster.indices.to.keep),] <- NA
        successes[!(clusters %in% cluster.indices.to.keep),] <- NA
        total.trials[!(clusters %in% cluster.indices.to.keep),] <- NA
      }

      r <- matrix(r[,indices.to.keep], nrow=N, ncol=N.c)
      ln.rho <- matrix(ln.rho[,indices.to.keep], nrow=N, ncol=N.c)
      # Need to renormalize r--do it gently.
      for(n in 1:N) {
         if(any(is.na(ln.rho[n,]))) {
           r[n,] <- rep(NA, N.c)
           next
         }

         row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
         for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
      }
      a <- matrix(a[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      b <- matrix(b[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      alpha <- alpha[indices.to.keep,drop=FALSE]
      a0 <- matrix(a0[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      b0 <- matrix(b0[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      alpha0 <- alpha0[indices.to.keep,drop=FALSE]

      E.pi.prev <- E.pi
    }

    if(do.inner.iteration == FALSE) { break }

  }

  retList <- list("retVal" = 0, "a" = a, "b" = b, "alpha" = alpha, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "outliers" = outliers)

  return(retList)

} # End binomial.bmm.filter.clusters


## ##--------------------------------------------------------------------------
## ## The binomial distribution clustering + filtering method
## ##
gaussian.bmm.filter.clusters <- function(vafs.merged, vafs, successes, total.trials, N.c, r, m, alpha, beta, nu, W, m0, alpha0, beta0, nu0, W0,
                                convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0,
                                plotIntermediateResults = 0, overlap.threshold=0.8){


  total.iterations <- 0
  num.dimensions <- dim(successes)[2]

  N <- dim(successes)[1]

  successes.colnames <- colnames(successes)
  total.colnames <- colnames(total.trials)

  outliers <- matrix(data=0, nrow=0, ncol=num.dimensions)
  colnames(outliers) <- colnames(vafs)

  E.pi.prev <- rep(0, N.c)

  width <- as.double(erf(1/sqrt(2)))

  while(TRUE) {

    if(verbose){
      print(r)
    }

    bmm.res <- gaussian.bmm.fixed.num.components(vafs, N.c, r, m, alpha, beta, nu, W, m0, alpha0, beta0, nu0, W0, convergence.threshold, max.iterations, verbose)
    if(bmm.res$retVal != 0) {
      cat("Failed to converge!\n")
      q(status=-1)
    }

    m <- bmm.res$m
    alpha <- bmm.res$alpha
    beta <- bmm.res$beta
    nu <- bmm.res$nu
    W <- bmm.res$W

    E.pi <- bmm.res$E.pi
    r <- bmm.res$r

    total.iterations <- total.iterations + bmm.res$num.iterations

    ln.rho <- bmm.res$ln.rho
    E.lnpi <- bmm.res$E.lnpi

    do.inner.iteration <- FALSE

    apply.min.items.condition <- TRUE
    apply.uncertainty.self.overlap.condition <- FALSE
    apply.large.SEM.condition <- FALSE
    apply.overlapping.SEM.condition <- FALSE
    apply.overlapping.std.dev.condition <- TRUE

    # Changes on Mar 25, 2013
    apply.overlapping.std.dev.condition <- FALSE
    apply.uncertainty.self.overlap.condition <- TRUE

    # remove.data = TRUE iff we should remove points assigned to clusters
    # that we remove
    remove.data <- FALSE

    clusters <- hardClusterAssignments(N,N.c,r)

    indices.to.keep <- rep(TRUE, N.c)
    if((apply.min.items.condition == TRUE) & (N.c > 1)) {
      threshold.pts <- 10

      num.items.per.cluster <- rep(0, N.c)
      for(n in 1:N) {
        num.items.per.cluster[clusters[n]] <- num.items.per.cluster[clusters[n]] + 1
      }

      indices.to.keep <- num.items.per.cluster >= threshold.pts

      # indices.to.keep <- E.pi > pi.threshold
    } # End apply.min.items.condition

    if(all(indices.to.keep==TRUE) & (apply.uncertainty.self.overlap.condition == TRUE) & (N.c > 1)) {

      overlaps <- rep(0, N.c)
      ones <- rep(1, N)
      # Just drop min overlap
      for(k in 1:N.c) {
        overlaps[k] <- sum(r[,k] * r[,k], na.rm=TRUE) / sum(ones * r[,k], na.rm=TRUE)
      }

      if(min(overlaps, na.rm=TRUE) < overlap.threshold) {
        for(k in 1:N.c) {
          # Small clusters will have undefined overlaps, just skip.
          # We'll remove them later.
          if(is.nan(overlaps[k])) { next }
          if((overlaps[k] < overlap.threshold) & (overlaps[k] == min(overlaps, na.rm=TRUE))) {
            indices.to.keep[k] <- FALSE
            break
          }
        }
      }

      if(verbose){
        for(k in 1:N.c) {
          cat(sprintf("Cluster %d pi = %.3f self-overlap = %.3f\n", k, E.pi[k], overlaps[k]))
        }

        means <- t(m)
        indices.to.drop <- (1:N.c)[!indices.to.keep]
        for(i in 1:length(indices.to.drop)) {
          index <- indices.to.drop[i]
          if(verbose){
            cat("Self-overlap condition: Dropping cluster with center: ")
            for(l in 1:num.dimensions) {
              cat(sprintf("%.3f ", means[l, index]))
            }
            cat("\n")
          }
        }
      }

    } # End apply.uncertainty.self.overlap.condition

    # Calculate standard error of the means
    SEM.res <- gaussian.bmm.narrowest.mean.interval.about.centers(m, alpha, beta, nu, W, width)
    SEM.centers <- t(SEM.res$centers)
    SEMs.lb <- t(SEM.res$lb)
    SEMs.ub <- t(SEM.res$ub)

    if(all(indices.to.keep==TRUE) & (apply.large.SEM.condition == TRUE) & (N.c > 1)) {

      for(k in 1:N.c) {
        if (verbose) {
          cat(sprintf("%dD: Cluster %d pi = %.3f: ", num.dimensions, k, E.pi[k]))
        }
        greater.than.30 <- TRUE
        greater.than.02 <- TRUE
        SEM.width.threshold <- 0.02
        #SEM.width.threshold <- 0.001
        for(m in 1:num.dimensions) {
          center <- SEM.centers[k,m]
          lower <- SEMs.lb[k,m]
          upper <- SEMs.ub[k,m]
          std.dev.width <- (std.dev.ub[k,m] - std.dev.lb[k,m])
          SEM.width <- (SEMs.ub[k,m] - SEMs.lb[k,m])
          #if((SEM.width/std.dev.width)>.3){ greater.than.30 <- TRUE }
          #if(SEM.width>.02){ greater.than.02 <- TRUE }
          if((SEM.width/std.dev.width)<=.3){ greater.than.30 <- FALSE }
          if(SEM.width<=SEM.width.threshold){ greater.than.02 <- FALSE }
          if (verbose) {
            cat(sprintf("%.3f (%.3f, %.3f) [(%.3f) %.3f, %.3f] {%.3f} ", center, lower, upper, width, std.dev.lb[k,m], std.dev.ub[k,m], SEM.width/std.dev.width))
            if(greater.than.30) { cat("* ") }
            if(greater.than.02) { cat("**") }
          }
        }
        if (verbose) {
          cat("\n")
        }
        if(greater.than.30 & greater.than.02) { indices.to.keep[k] <- FALSE }
      }
      if ( any(indices.to.keep==FALSE) ) {
        remove.data <- TRUE
      }
    } # End apply.large.SEM.condition

    if(all(indices.to.keep==TRUE) & (apply.overlapping.SEM.condition == TRUE) & (N.c > 1)) {

      # Determine if component i's center is contained within
      # component i2's std.dev
      pi.threshold = 10^-2

      for(i in 1:N.c){
        i.subsumed.by.another.cluster <- FALSE
        for(i2 in 1:N.c){
          if(i == i2) { next }
          if(indices.to.keep[i2] == FALSE) { next }
          if(E.pi[i2] < pi.threshold) { next }
          i.subsumed.by.another.cluster <- TRUE
          for(l in 1:num.dimensions){
            i.center <- std.dev.centers[i,l]
            if((i.center < std.dev.lb[i2,l]) | (i.center > std.dev.ub[i2,l])) {
              i.subsumed.by.another.cluster <- FALSE
              break
            }
          }
          if(i.subsumed.by.another.cluster==TRUE) {
            if(verbose){
              cat(sprintf("2. Dropping cluster with center: "))
              for(l in 1:num.dimensions){
                cat(sprintf("%.3f ", std.dev.centers[i,l]))
              }
              cat(sprintf("because it overlaps with: "))
              for(l in 1:num.dimensions){
                cat(sprintf("(%.3f, %.3f) ", std.dev.lb[i2,l], std.dev.ub[i2,l]))
              }
              cat("\n")
              break
            }
          }
        }
        if(i.subsumed.by.another.cluster==TRUE) {
          indices.to.keep[i] <- FALSE
        }
      }

    } # End apply.overlapping.SEM.condition

    if(all(indices.to.keep==TRUE) & (apply.overlapping.std.dev.condition == TRUE) & (N.c > 1)) {

      # Determine how much component i's std dev overlaps component i2's std dev
      overlaps <- matrix(data = 1, nrow=N.c, ncol=N.c)
      std.dev.overlap.threshold <- 0
      for(i in 1:N.c){
        for(i2 in 1:N.c){
          if(verbose){
            cat("   ")
          }
          for(l in 1:num.dimensions){
            fraction <- 0
            if((std.dev.lb[i,l] < std.dev.ub[i2,l]) & (std.dev.ub[i,l] > std.dev.lb[i2,l])) {
              overlap <- min(std.dev.ub[i,l], std.dev.ub[i2,l]) - max(std.dev.lb[i,l], std.dev.lb[i2,l])
              fraction <- overlap / ( std.dev.ub[i,l] - std.dev.lb[i,l] )
            }
            if((fraction > std.dev.overlap.threshold) & (overlaps[i,i2] != 0)) {
              overlaps[i,i2] <- 1
            } else {
              overlaps[i,i2] <- 0
            }
          }
        }
      }

      # Remove one of two overlapping clusters.  If both overlap
      # (NB: above is not symmetric), then only remove smaller of two.
      for(i in 2:N.c){
        if(indices.to.keep[i] == FALSE) { next }
        for(i2 in 1:(i-1)){
          if(indices.to.keep[i2] == FALSE) { next }
          if(overlaps[i,i2] == 1) {
            if((overlaps[i2,i] == 1) & (E.pi[i2] < E.pi[i])) {
              if(verbose){cat("Removing ", i2, " because of overlap with i = ", i, "\n")}
              indices.to.keep[i2] <- FALSE
            } else {
              if(verbose){cat("Removing ", i, " because of overlap with i2 = ", i2, "\n")}
              indices.to.keep[i] <- FALSE
            }
          }
        }
      }
    } # End apply.overlapping.std.dev.condition



    if ( any(indices.to.keep==FALSE) ) {

      do.inner.iteration <- TRUE

      numeric.indices <- (1:N.c)

      E.pi <- E.pi[indices.to.keep]
      E.lnpi <- E.lnpi[indices.to.keep]

      N.c <- length(E.pi)

      if(remove.data == TRUE) {
        cluster.indices.to.keep <- (1:N.c)[indices.to.keep]
        new.outliers <- matrix(vafs[!(clusters %in% cluster.indices.to.keep),], ncol=num.dimensions)
        outliers <- rbind(outliers, new.outliers)
        # To remove data from the data set, set its entries to NA
        vafs[!(clusters %in% cluster.indices.to.keep),] <- NA
        successes[!(clusters %in% cluster.indices.to.keep),] <- NA
        total.trials[!(clusters %in% cluster.indices.to.keep),] <- NA
      }

      r <- matrix(r[,indices.to.keep], nrow=N, ncol=N.c)
      ln.rho <- matrix(ln.rho[,indices.to.keep], nrow=N, ncol=N.c)
      # Need to renormalize r--do it gently.
      for(n in 1:N) {
         if(any(is.na(ln.rho[n,]))) {
           r[n,] <- rep(NA, N.c)
           next
         }

         row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
         for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
      }

      m <- matrix(m[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      alpha <- alpha[indices.to.keep,drop=FALSE]
      beta <- beta[indices.to.keep,drop=FALSE]
      nu <- nu[indices.to.keep,drop=FALSE]
      W <- W[indices.to.keep,drop=FALSE]

      m0 <- matrix(m0[indices.to.keep,], nrow=N.c, ncol=num.dimensions)
      alpha0 <- alpha0[indices.to.keep,drop=FALSE]
      beta0 <- beta0[indices.to.keep,drop=FALSE]
      nu0 <- nu0[indices.to.keep,drop=FALSE]
      W0 <- W0[indices.to.keep,drop=FALSE]

      E.pi.prev <- E.pi
    }

    if(do.inner.iteration == FALSE) { break }

  }

  L <- gaussian.bmm.calculate.posterior.predictive.precision(m, alpha, beta, nu, W)

  retList <- list("retVal" = 0, "m" = m, "alpha" = alpha, "beta" = beta, "nu" = nu, "W" = W, "L" = L, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "outliers" = outliers)

  return(retList)

} # End gaussian.bmm.filter.clusters
