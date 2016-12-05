# Based on adaptive SMC algorithm proposed by Del Moral, Pierre, Arnaud
# Doucet, and Ajay Jasra. "An adaptive sequential Monte Carlo method for
# approximate Bayesian computation." Statistics and Computing 22.5 (2012):
# 1009-1020.


# global parameters
resize.amount <- 100
bisection.max.iter <- 10000


simulate.trees <- function(workspace, theta, seed=NA, ...) {
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
    result <- config$model(theta, config$nsample, workspace$n.tips,
        workspace$tip.heights, workspace$tip.labels, ...)

    # annotate each trees with its self-kernel score
    for (i in 1:config$nsample) {
        result[[i]]$kernel <- tree.kernel(
            result[[i]],
            result[[i]],
            lambda=config$decay.factor,
            sigma=config$rbf.variance,
            rho=config$sst.control,
            rescale.mode=config$norm.mode
        )
    }

    return(result)
}


distance <- function(t1, t2, config) {
    k <- tree.kernel(
        sim.tree,
        obs.tree,
        lambda=config$decay.factor,
        sigma=config$rbf.variance,
        rho=config$sst.control,
        rescale.mode=config$norm.mode
    )
    return (1. - k / sqrt(sim.tree$kernel * obs.tree$kernel))
}


initialize.smc <- function(ws, ...) {
    config <- ws$config
    nparams <- len(config$params)

	for (i in 1:config@nparticle) {
        # sample particle from prior distribution
		ws$particles[i,] <- sample(config)

        # assign uniform weights
		ws$weights[i] <- 1./config@nparticle

		# simulate trees from particle
		ws$sim.trees[[i]] <- simulate.trees(ws, config, ...)

		# calculate kernel distances for trees
		ws$dists[,i] <- sapply(ws$sim.trees[[i]], function(sim.tree) {
            distance(obs.tree, sim.tree, config)
		})
	}
    cat('Initialized SMC workspace.\n')
}



ess <- function(w) {
    # effective sample size
    return(1/sum(w^2))
}

epsilon.obj.func <- function(ws, epsilon) {
    # unpack some things
    config <- ws$config
    prev.epsilon <- ws$epsilon

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
    return (ess(ws$new.weights) - config$alpha * ess(ws$weights))
}


next.epsilon <- function(ws) {
    # Let W_n^i be the weight of the i-th particle at n-th iteration
    #
    # The effective sample size is
    #   ESS({W_n^i}) = 1 / \sum_{i=1}^{N} (W_n^i)^2
	# Use bisection method to solve for epsilon such that:
    #   ESS(W*, eps) - alpha * ESS(W, eps) = 0
    config <- ws$config
    res <- uniroot(function(x) epsilon.obj.func(ws, x), lower=0,
        upper=ws$epsilon, tol=config$step.tolerance)
    root <- res$root
    if (root < config$final.epsilon) {
        root = config$final.epsilon  # stopping criterion
        epsilon.obj.func(ws, root)
    }
    ws$weights <- ws$new.weights  # update weights
    return (root)
}


resample.particles <- function(ws) {
}


perturb.particles <- function(ws) {
}

run.smc <- function(ws, trace.file=NA, regex=NA, seed=NA, nthreads=1, ...) {
    # @param ws: workspace
	# @param obs.tree: object of class 'phylo'
	# @param trace.file: (optional) path to a file to write outputs
	# @param seed: (optional) integer to set random seed
	# @param nthreads: (optional) for running on multiple cores
	# @param ...: additional arguments to pass to config@generator via
    #   simulate.trees()

	config <- ws$config

    # space for returned values
	accept.rate <- {}
	epsilons <- {}
	
    initialize.smc(ws)

    n.iter <- 0
    epsilon <- .Machine$double.xmax
    while (epsilon != config$final.epsilon) {
        ws$accept <- 0
        ws$alive <- 0

        # update epsilon
        ws$epsilon <- next.epsilon(ws)

        # resample particles according to their weights
        if (ess(ws$weights) < config$ess.tolerance) {
            resample.particles(ws)
        }

        # perturb particles
        ws$accept <- 0
        ws$alive <- 0
        perturb.particles(ws)
    }
}



