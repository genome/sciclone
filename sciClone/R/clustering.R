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
##   fit.x = a vector of length P, where P is an arbitrary number of points between 0 and 1
##           where the model fit was sampled (in each of the M dimensions)
##   fit.y = a matrix of size MxP, containing the corresponding Y value for each X
##           Y values should be scaled between 0 and 1
##   individual.fits.y = a list of length number_of_clusters, holding the individual fits for each of the models, each represented by an M x P matrix (as for fit.y)

clusterVafs <- function(vafs.merged, vafMatrix, method="bmm", purities=100, params=NULL, samples=1, plotIntermediateResults = 0, verbose=0){
  ##check for suitable method
  if(method == "bmm"){
   return(clusterWithBmm(vafs.merged, vafMatrix, samples=samples, plotIntermediateResults=plotIntermediateResults, verbose=0))
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
clusterWithBmm <- function(vafs.merged, vafs, initialClusters=10, samples=1, plotIntermediateResults=0, verbose=TRUE) {
    library(bmm)

    #replace any values of zero with a very small number to prevent errors
    delta <- .Machine$double.eps
    vafs[which(vafs==0)] = delta

    ## Initialize the hyperparameters of the Beta mixture model (bmm).
    ## print(length(vafs))
    ## print(head(vafs))
    ## print(head(vafs.merged))
    hyperparams <- init.bmm.hyperparameters(vafs, initialClusters)

    ## Initialize the parameters of the bmm.
    params <- init.bmm.parameters(vafs, initialClusters, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0)

    ## Perform the clustering.
    ## Start with the provided number of clusters, but prune any with low probability
    bmm.results <- bmm.filter.clusters(vafs.merged, vafs, initialClusters, params$r, params$mu, params$alpha, params$nu, params$beta, params$c, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose, plotIntermediateResults=plotIntermediateResults)
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


## ##--------------------------------------------------------------------------
## ## The beta distribution clustering + filtering method
## ##
bmm.filter.clusters <- function(vafs.merged, X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4,
                                max.iterations = 10000, verbose = 0, plotIntermediateResults = 0){

    total.iterations <- 0
    N <- dim(X)[1] 
    num.dimensions <- dim(X)[2]

    x.colnames <- colnames(X) 

    outliers <- matrix(data=0, nrow=0, ncol=dim(X)[2])
    colnames(outliers) <- colnames(X)

    E.pi.prev <- rep(0, N.c)

    width <- as.real(erf(1.5/sqrt(2)))
    # width <- as.real(erf(1/sqrt(2)))
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

      ellipse.width <- as.real(erf(1/sqrt(2)))
      
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

          ellipse.width <- as.real(erf(1/sqrt(2)))
      
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
        apply.uncertainty.overlap.condition <- TRUE        

        if((apply.min.items.condition == TRUE) & (N.c > 1)) {
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

            overlap.threshold <- 0.8
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
