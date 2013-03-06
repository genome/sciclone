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
## Go from fuzzy probabilities to hard cluster assignments
##
hardClusterAssignments <- function(numPoints,numClusters,probabilities) {

    print("Inside Hard Cluster Assignments");
    print(numPoints);
    print(numClusters);
    print(dim(probabilities));

    assignments <- rep(NA,numPoints);
    for(n in 1:numPoints) {
        max.cluster <- 0
        max.assignment <- -1
        for(k in 1:numClusters) {
            if ( probabilities[n,k] > max.assignment ) {
                max.assignment <- probabilities[n,k]
                max.cluster <- k
            }
        }
        assignments[n] <- max.cluster
    }
    return(assignments)
}


##--------------------------------------------------------------------------
## Do clustering with bmm (binomial mixture model)
##
clusterWithBmm <- function(vafs, initialClusters=10, samples=1) {
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
    bmm.results <- bmm.filter.clusters(vafs, initialClusters, params$r, params$mu, params$alpha, params$nu, params$beta, params$c, hyperparams$mu0, hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, hyperparams$c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
    if(bmm.results$retVal != 0) {
        cat("WARNING: bmm failed to converge. No clusters assigned\n")
        return(NULL);
    }

    ##get the assignment of each point to a cluster
    probs = exp(1)^bmm.results$ln.rho
    numPoints = length(probs[,1])
    numClusters = length(probs[1,])
    clusters = hardClusterAssignments(numPoints,numClusters,probs);

    ## find confidence intervals around the means of the clusters
    intervals = bmm.narrowest.mean.interval.about.centers(bmm.results$mu, bmm.results$alpha, bmm.results$nu, bmm.results$beta, 0.68)
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
                ym[k,i] <- bmm.component.posterior.predictive.density(x[i], bmm.results$mu[dim,k], bmm.results$alpha[dim,k], bmm.results$nu[dim,k], bmm.results$beta[dim,k], bmm.results$E.pi[k], num.samples = num.iterations)
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
## ## The beta distribution clustering + filtering method
## ##
bmm.filter.clusters <- function(X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
{

    total.iterations <- 0
    N <- dim(X)[1] 
    num.dimensions <- dim(X)[2]

    x.colnames <- colnames(X) 

    outliers <- matrix(data=0, nrow=0, ncol=dim(X)[2])
    colnames(outliers) <- colnames(X)

    E.pi.prev <- rep(0, N.c)

    width <- as.real(erf(1.5/sqrt(2)))

    # Outer while loop: following convergence of inner loop, apply
    # overlapping cluster condition to drop any overlapping clusters.
    while(TRUE) {

        print(dim(X))
        print(N.c)
        print(N)
        print(length(c0))
        print(dim(r))
        print(length(colSums(r)))

        bmm.res <- bmm.fixed.num.components(X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0, convergence.threshold, max.iterations, verbose)
        if(bmm.res$retVal != 0) {
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

        print(dim(r))
        print(length(colSums(r)))

        total.iterations <- total.iterations + bmm.res$num.iterations

        ln.rho <- bmm.res$ln.rho
        E.lnu <- bmm.res$E.lnu
        E.lnv <- bmm.res$E.lnv
        E.lnpi <- bmm.res$E.lnpi
        E.quadratic.u <- bmm.res$E.quadratic.u
        E.quadratic.v <- bmm.res$E.quadratic.v

        do.inner.iteration <- FALSE

        # Remove any small clusters

        apply.min.items.condition <- TRUE
        apply.uncertainty.condition <- FALSE
        apply.large.SEM.condition <- FALSE
        apply.overlapping.SEM.condition <- FALSE
        apply.overlapping.std.dev.condition <- TRUE

        print("inside bmm, before apply.min.items.cond if statement");
        print(paste("apply.min.items = ",apply.min.items.condition));
        print(paste("N.c =",N.c));


        if((apply.min.items.condition == TRUE) & (N.c > 1)) {
            threshold.pts <- 3

            print("inside bmm, AFTER apply.min.items.cond if statement");
        print(paste("apply.min.items = ",apply.min.items.condition));
        print(paste("N.c =",N.c));
            
            clusters <- hardClusterAssignments(N,N.c,r)

            num.items.per.cluster <- rep(0, N.c)
            for(n in 1:N) {
                num.items.per.cluster[clusters[n]] <- num.items.per.cluster[clusters[n]] + 1
            }

            non.zero.indices <- num.items.per.cluster >= threshold.pts

            if ( any(non.zero.indices==FALSE) ) {

                do.inner.iteration <- TRUE

                numeric.indices <- (1:N.c)
                cat("Dropping clusters with pts: \n")
                print(num.items.per.cluster[!non.zero.indices])  

                remove.data <- TRUE
                if (remove.data == TRUE) {
                    numeric.indices <- (1:N.c)[non.zero.indices]
                    new.outliers <- matrix(X[!(clusters %in% numeric.indices),], ncol=num.dimensions)
                    outliers <- rbind(outliers, new.outliers)
                }

                X <- matrix(X[clusters %in% numeric.indices,], ncol=num.dimensions)
                colnames(X) <- x.colnames
                N <- dim(X)[1]

                E.pi <- E.pi[non.zero.indices]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[non.zero.indices]
                c <- c[non.zero.indices]
                c0 <- c0[non.zero.indices]

                r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
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

        if((apply.uncertainty.condition == TRUE) & (N.c > 1)) {

            overlaps <- rep(0, N.c)
            ones <- rep(1, N)
            indices.to.keep.boolean <- rep(TRUE, N.c)
            # Just drop min overlap
            for(k in 1:N.c) {
                overlaps[k] <- sum(r[,k] * r[,k]) / sum(ones * r[,k])
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

            for(k in 1:N.c) {
                # cat(sprintf("Cluster %d pi = %.3f self-overlap = %.3f\n", num.dimensions, k, E.pi[k], overlaps[k]))
            }

            indices.to.keep <- (1:N.c)[indices.to.keep.boolean]

            if(length(indices.to.keep) != N.c) {

                do.inner.iteration <- TRUE

                means <- matrix(ubar/(ubar+vbar), nrow=num.dimensions, ncol=N.c)
                if(length(indices.to.keep) < N.c) {
                    indices.to.drop <- (1:N.c)[!indices.to.keep.boolean]
                    for(i in 1:length(indices.to.drop)) {
                        index <- indices.to.drop[i] 
                        cat("1. Dropping cluster with center: ")
                        for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                        }
                        cat("\n")
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
        } # End apply.uncertainty.condition

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

            show.output <- TRUE
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

                cat("Removing clusters because of SEM condition: \n")
                print(numeric.indices[non.zero.indices])
                #HERE TEST
                clusters <- rep(0, N)
                for(n in 1:N) {
                    max.cluster <- 0
                    max.assignment <- -1
                    for(k in 1:N.c) {
                        if ( r[n,k] > max.assignment ) {
                            max.assignment <- r[n,k]
                            max.cluster <- k
                        }
                    }
                    clusters[n] <- max.cluster
                }

                remove.data <- TRUE
                if (remove.data == TRUE) {
                    numeric.indices <- (1:N.c)[non.zero.indices]
                    new.outliers <- matrix(X[!(clusters %in% numeric.indices),], ncol=num.dimensions)
                    outliers <- rbind(outliers, new.outliers)
                }

                X <- matrix(X[clusters %in% numeric.indices,], ncol=num.dimensions)
                colnames(X) <- x.colnames
                N <- dim(X)[1]

                E.pi <- E.pi[non.zero.indices]
                N.c <- length(E.pi)
                E.pi.prev <- E.pi.prev[non.zero.indices]
                c <- c[non.zero.indices]
                c0 <- c0[non.zero.indices]

                r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
                ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
                # Need to renormalize r--do it gently.
                for(n in 1:N) {
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
                        cat(sprintf("2. Dropping cluster with center: "))
                        for(l in 1:num.dimensions){
                            cat(sprintf("%.3f ", std.dev.centers[i,l]))
                        }
                        cat(sprintf("becuase it overlaps with: "))
                        for(l in 1:num.dimensions){
                            cat(sprintf("(%.3f, %.3f) ", std.dev.lb[i2,l], std.dev.ub[i2,l]))
                        }
                        cat("\n")
                        break
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
                        cat("3. Dropping cluster with center: ")
                        for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                        }
                        cat("\n")
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
                    cat("   ")
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
                            cat("Removing ", i2, " because of overlap with i = ", i, "\n")
                            indices.to.keep.boolean[i2] <- FALSE               
                        } else {
                            cat("Removing ", i, " because of overlap with i2 = ", i2, "\n")
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
                        cat("4. Dropping cluster with center: ")
                        for(l in 1:num.dimensions) {
                            cat(sprintf("%.3f ", means[l, index]))
                        }
                        cat("\n")
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
