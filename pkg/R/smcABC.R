# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.


# Based on adaptive SMC algorithm proposed by Del Moral, Pierre, Arnaud
# Doucet, and Ajay Jasra. "An adaptive sequential Monte Carlo method for
# approximate Bayesian computation." Statistics and Computing 22.5 (2012):
# 1009-1020.


# global parameters
resize.amount <- 100
bisection.max.iter <- 10000


simulate.trees <- function(workspace, theta, model, seed=NA, nthreads=1, ...) {
	# @param workspace: smc.workspace object
	# @param theta: parameter vector
	# @param seed: argument to set.seed()
  config <- workspace$config
  if (is.null(body(config$model))) {
    stop('Simulation method has not been set for configuration.')
  }
	if (!is.na(seed)) {
		set.seed(seed)
	}
  result <- config$model(theta=theta, nsim=config$nsample, tips=workspace$n.tips,
                         model=model, seed=seed, labels=workspace$tip.labels, ...)
  
  # annotate each tree with its self-kernel score (parallelized)
  processed <- mclapply(1:config$nsample,
                        function(i) {
                        result[[i]] <<- .preprocess.tree(result[[i]], config)
                        },
                        mc.cores = nthreads)
  
  return(processed)
}


# formally 'distance'
kernel.dist <- function(t1, t2, decay.factor, rbf.variance, sst.control, norm.mode) {
  if (is.null(t1$kernel)) {
    stop("t1 missing self kernel in distance()")
  }
  if (is.null(t2$kernel)) {
    stop("t2 missing self kernel in distance()")
  }
  
  k <- tree.kernel(
    t1,
    t2,
    lambda=decay.factor,
    sigma=rbf.variance,
    rho=sst.control
  )
  
  result <- 1. - k / sqrt(t1$kernel * t2$kernel)
  if (result < 0 || result > 1) {
    stop(
      cat("ERROR: kernel.dist() value outside range [0,1].\n",
          "k: ", k, "\n",
          "t1$kernel: ", t1$kernel, "\n",
          "t2$kernel: ", t2$kernel, "\n"
      )
    )
  }
  if (is.nan(result)) {
    cat("t1$kernel:", t1$kernel, "\n")
    cat("t2$kernel:", t2$kernel, "\n")
  }
  return (result)
}


# Applies config$dist expression to trees x and y
distance <- function(x, y, config) {
  
  result <- eval(parse(text=config$dist)) 
  
  # One instance of NA or NAN throws error
  if (is.na(result)==TRUE) {
    stop("One or more distance metrics returns NA")
  } else if (is.nan(result)==TRUE) {
    stop("One or more distance metrics returns NAN")
  } else {
    return(result)
  }
}


initialize.smc <- function(ws, model, seed=NA, nthreads=1, ...) {
  config <- ws$config
  nparams <- length(config$params)
  colnames(ws$particles) <- config$params
  for (i in 1:config$nparticle) {
    # sample particle from prior distribution
    ws$particles[i,] <- sample.priors(config)
    
    # assign uniform weights
    ws$weights[i] <- 1./config$nparticle
    
    # simulate trees from particle
    ws$sim.trees[[i]] <- simulate.trees(ws, ws$particles[i,], model=model, seed=seed, nthreads=nthreads, ...)
    
		# calculate kernel distances for trees (parallelized)
		ws$dists[,i] <- unlist(mclapply(ws$sim.trees[[i]], 
		                                function(sim.tree) {
		                                distance(ws$obs.tree, sim.tree, config)
		                                },
		                                mc.cores=nthreads))
	}
  cat('Initialized SMC workspace.\n')
  return(ws)
}


.ess <- function(w) {
  # effective sample size
  return(1/sum(w^2))
}


.epsilon.obj.func <- function(ws, epsilon) {
  # unpack some things
  config <- ws$config
  prev.epsilon <- ws$epsilon

  # check that all required objects are present
  if (is.null(config$quality)) {
    stop("epsilon.obj.func: Config missing setting `quality`, exiting")
  }

  # calculate new weights
  for (i in 1:config$nparticle) {
    num <- sum(ws$dists[,i] < epsilon)
    denom <- sum(ws$dists[,i] < prev.epsilon)
    
    if (num == denom) {
      # handle case where numerator and denominator are both zero
      ws$new.weights[i] <- ws$weights[i]
    } else {
      ws$new.weights[i] <- ws$weights[i] * num / denom
    }
  }

  # normalize new weights to sum to 1.
  wsum <- sum(ws$new.weights)
  ws$new.weights <- ws$new.weights / wsum
  if (epsilon==0 || wsum==0) {
    return (-1)  # undefined -- range limited to (0, prev.epsilon)
  }
  return (.ess(ws$new.weights) - config$quality * .ess(ws$weights))
}


.next.epsilon <- function(ws, niter) {
  # Let W_n^i be the weight of the i-th particle at n-th iteration
  #
  # The effective sample size is
  #   ESS({W_n^i}) = 1 / \sum_{i=1}^{N} (W_n^i)^2
	# Use bisection method to solve for epsilon such that:
  #   ESS(W*, eps) - alpha * ESS(W, eps) = 0, where alpha is quality parameter

  # check that dists have been set
  if (any(is.na(ws$dists))) {
    stop("NA values in ws$dists; did you forget to run initialize.smc()?")
  }
  config <- ws$config
  # solve for new epsilon
  res <- uniroot(function(x) .epsilon.obj.func(ws, x), lower=0,
                 upper=ws$epsilon, tol=config$step.tolerance, maxiter=1E6)
  root <- res$root
  if (root < config$final.epsilon) {
    cat("adjusted root from ", root, " to ", config$final.epsilon, "\n")
    root = config$final.epsilon
  }

  # recalculate new weights at new epsilon
  #   this is annoying but it's difficult in R to pass `ws` by reference
  #   to epsilon.obj.func where this has been done already...
  for (i in 1:config$nparticle) {
    num <- sum(ws$dists[,i] < root)
    denom <- sum(ws$dists[,i] < ws$epsilon)
    ws$weights[i] <- ws$weights[i] * ifelse(num==denom, 1., num/denom)
    
  }
  ws$weights <- ws$weights / sum(ws$weights)  # renormalize weights
  ws$epsilon <- root
  return (ws)
}


.resample.particles <- function(ws) {
  nparticle <- ws$config$nparticle
  # sample from current population of particles with replacement
  indices <- sample(1:nparticle, nparticle, replace=TRUE, prob=ws$weights)

  # if there's only one parameter, this returns a vector unless we recast
  ws$particles <- as.matrix(ws$particles[indices,])
  colnames(ws$particles) <- ws$config$params

  ws$dists <- ws$dists[,indices]  # transfer columns of kernel distances

  # reset all weights
  ws$weights <- rep(1./nparticle, times=nparticle)
  return(ws)
}



.perturb.particles <- function(ws, model, nthreads=1, seed=NA, ...) {
  ##  This implements the Metropolis-Hastings acceptance/rejection step
  ##  @param ws: workspace
  ##  @param model: reference to an R function that will simulate trees
  ##  @param nthreads:  number of threads to run in parallel, defaults to 1 (serial)
  if (nthreads < 1) {
    stop("User requested less than 1 thread in call to .perturb.particles")
  }

  config <- ws$config
  nparticle <- config$nparticle

  require(parallel, quietly=TRUE)  # load parallel library if not already present

  # iterate over live particles
  alive <- which(ws$weights > 0)  # indices of live particles
  ws$alive <- length(alive)
  
  res <- mclapply(alive, function (i) {
    old.particle <- ws$particles[i,]
    new.particle <- propose(config, old.particle)
    
    # calculate prior ratio
    mh.ratio <- prior.density(config, new.particle) / prior.density(config, old.particle)
    if (mh.ratio == 0) {
      return(NULL)    # reject new particle, violates prior assumptions
    }
    
    # calculate proposal ratio
    mh.ratio <- mh.ratio * proposal.density(config, new.particle, old.particle) / proposal.density(config, old.particle, new.particle)
    if (mh.ratio == 0) {
      return(NULL)    # reject new particle, not possible under proposal distribution
    }
    
    # simulate new trees (parallelized)
    # retain sim.trees in case we revert to previous particle
    new.trees <- simulate.trees(ws, new.particle, model=model, seed=seed, nthreads=nthreads, ...)
    new.dists <- unlist(mclapply(new.trees, 
                                    function(sim.tree) {
                                      distance(ws$obs.tree, sim.tree, config)
                                    },
                                    mc.cores=nthreads))
    
    # SMC approximation to likelihood ratio
    old.nbhd <- sum(ws$dists[,i] < ws$epsilon)    # how many samples are in the neighbourhood of data?
    new.nbhd <- sum(new.dists < ws$epsilon)
    mh.ratio <- mh.ratio * new.nbhd / old.nbhd
    
    output <- list(i, mh.ratio, new.particle, new.trees, new.dists)
    output
    
  }, mc.cores=nthreads)  

  # accept or reject the proposal
  for (row in res) {
    if (length(row) == 0) {            # checking for any items that are a returned NULL
      next
    }
    else {
      i <- row[[1]]
      if (runif(1) < row[[2]]) {     # always accept if mh.ratio > 1 
        
        ws$accept[i] <- TRUE
        ws$particles[i,] <- row[[3]]    # new.particle
        ws$dists[,i] <- row[[5]]        # new.dists             
        ws$sim.trees[[i]] <- row[[4]]   # new.trees
      }
      else {
        ws$accept[i] <- FALSE
      }
    }
  }
  # creating new attribute ws$accepted in workspace; didn't want dual vector and int behaviour of ws$accept from parallelization
  ws$accepted <- length(which(ws$accept == TRUE))  # FIXME: <- sum(ws$accept)  # should work..
  # TODO: use return values to update ws
  return (ws)
}



run.smc <- function(ws, trace.file='', regex=NA, seed=NA, nthreads=1, verbose=FALSE, model='', maxReject=10, ...) {
  # @param ws: workspace
	# @param obs.tree: object of class 'phylo'
	# @param trace.file: (optional) path to a file to write outputs
	# @param seed: (optional) integer to set random seed
	# @param nthreads: (optional) for running on multiple cores
	# @param ...: additional arguments to pass to config@generator via
  #   simulate.trees()
  
  config <- ws$config
  
  # clear file and write header row
  write.table(t(c(
    'n', 'part.num', 'weight', config$params, paste0('dist.', 1:config$nsample)
    )), file=trace.file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

  # space for returned values
  result <- list(niter=0, theta=list(), weights=list(), accept.rate={}, epsilons={})

  # draw particles from prior distribution, assign weights and simulate data
  ptm <- proc.time()  # start timer
  cat ("Initializing SMC-ABC run with", config$nparticle, "particles\n")
  ws <- initialize.smc(ws, model, seed=seed, nthreads=nthreads, ...)

  niter <- 0
  ws$epsilon <- .Machine$double.xmax

  # report stopping conditions
  if (verbose) {
    cat ("ws$epsilon: ", ws$epsilon, "\n");
    cat ("config$final.epsilon: ", config$final.epsilon, "\n");
  }

  while (ws$epsilon != config$final.epsilon) {
    niter <- niter + 1

    # update epsilon
    ws <- .next.epsilon(ws, niter)

    # provide some feedback
    lap <- proc.time() - ptm
    cat ("Step ", niter, " epsilon:", ws$epsilon, " ESS:", .ess(ws$weights),
             "accept:", result$accept.rate[length(result$accept.rate)],
             "elapsed:", round(lap[['elapsed']],1), "s\n")

    # resample particles according to their weights
    if (.ess(ws$weights) < config$ess.tolerance) {
      ws <- .resample.particles(ws)
    }

    # perturb particles
    ws$accept <- vector()    # vector to keep track of which particles were accepted through parallelization in .perturb.particles
    ws$alive <- 0
    ws <- .perturb.particles(ws, model, nthreads=nthreads)  # Metropolis-Hastings sampling
    
    # record everything
    result$theta[[niter]] <- ws$particles
    result$weights[[niter]] <- ws$weights
    result$epsilons <- c(result$epsilons, ws$epsilon)
    result$accept.rate <- c(result$accept.rate, ws$accepted / ws$alive)      # changed ws$accept to ws$accepted; didn't want dual behaviour of ws$accept switching back and forth between vector and int

    # write output to file if specified
    for (i in 1:config$nparticle) {
      write.table(
                x=t(c(niter, i, round(ws$weights[i],10), round(ws$particles[i,],5), round(ws$dists[,i], 5))),
                file=trace.file,
                append=TRUE,
                sep="\t",
                row.names=FALSE,
                col.names=FALSE
            )
    }

    # report stopping conditions
    if (verbose) {
      cat("run.smc niter: ", niter, "\n")
      cat ("ws$epsilon: ", ws$epsilon, "\n");
      cat ("config$final.epsilon: ", config$final.epsilon, "\n");
      cat ("result$accept.rate: ", result$accept.rate, "\n");
      cat ("config$final.accept.rate: ", config$final.accept.rate, "\n");
    }

    # if acceptance rate is low enough, we're done
    if (result$accept.rate[niter] <= config$final.accept.rate) {
      ws$epsilon <- config$final.epsilon
      break  # FIXME: this should be redundant given loop condition above
    }
    
    ## check if the last 10 accept.rates OR last 10 epislon values are the same: if frozen, break
    if (niter > maxReject &&
        (length(unique(result$accept.rate[niter:(niter-maxReject-1)])) == 1 
         || length(unique(result$epsilons[niter:(niter-maxReject-1)])) == 1)) {
        cat ("SMC-ABC run has frozen in given parameter space. Please check the ranges in your prior settings.\n")
        break
    }
  }

  # finally sample from the estimated posterior distribution
  ws <- .resample.particles(ws)
  result$theta[[niter]] <- ws$particles
  result$weights[[niter]] <- ws$weights
  result$niter <- niter
  
  for (i in 1:config$nparticle) {
    write.table(
      x=t(c((niter + 1), i, round(ws$weights[i],10), round(ws$particles[i,],5), round(ws$dists[,i], 5))),
      file=trace.file,
      append=TRUE,
      sep="\t",
      row.names=FALSE,
      col.names=FALSE
    )
  }  
  # pack ws and result into one list to be returned
  ret <- list(workspace=ws, result=result)
    
  return (ret)
}



