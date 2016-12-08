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
        t1,
        t2,
        lambda=config$decay.factor,
        sigma=config$rbf.variance,
        rho=config$sst.control,
        rescale.mode=config$norm.mode
    )
    return (1. - k / sqrt(t1$kernel * t2$kernel))
}


initialize.smc <- function(ws, ...) {
    config <- ws$config
    nparams <- length(config$params)
    colnames(ws$particles) <- config$params

	for (i in 1:config$nparticle) {
        # sample particle from prior distribution
		ws$particles[i,] <- sample.priors(config)

        # assign uniform weights
		ws$weights[i] <- 1./config$nparticle

		# simulate trees from particle
		ws$sim.trees[[i]] <- simulate.trees(ws, ws$particles[i,], ...)

		# calculate kernel distances for trees
		ws$dists[,i] <- sapply(ws$sim.trees[[i]], function(sim.tree) {
            distance(ws$obs.tree, sim.tree, config)
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
    nparticle <- ws$config$nparticle
    # sample from current population of particles with replacement
    indices <- sample(1:nparticle, nparticle, replace=TRUE, prob=ws$weights)
    ws$particles <- ws$particles[,indices]
    ws$dists <- ws$dists[,indices]  # transfer columns of kernel distances

    # reset all weights
    ws$weights <- rep(1./nparticle, times=nparticle)
    return(ws)
}


perturb.particles <- function(ws) {
    ##  This implements the Metropolis-Hastings acceptance/rejection step
    config <- ws$config
    nparticle <- config$nparticle

    # loop over particles
    # TODO: multi-threaded implementation
    for (i in 1:nparticle) {
        if (ws$weights[i] == 0) {
            next  # ignore dead particles
        }
        ws$alive <- ws$alive + 1
        old.particle <- ws$particle[i,]
        new.particle <- propose(config, old.particle)

        # calculate prior ratio
        mh.ratio <- prior.density(config, new.particle) / prior.density(config, old.particle)
        if (mh.ratio == 0) {
            next  # reject new particle, violates prior assumptions
        }

        # calculate proposal ratio
        mh.ratio <- mh.ratio * proposal.density(config, new.particle, old.particle) /
                    proposal.density(config, old.particle, new.particle)
        if (mh.ratio == 0) {
            next  # reject new particle, not possible under proposal distribution
        }

        # simulate new trees  # TODO: this is probably a good spot for parallelization
        for (j in 1:config$nsample) {
            # retain sim.trees in case we revert to previous particle
            new.trees <- simulate.trees(ws, new.particle)
            new.dists <- sapply(new.trees, function(sim.tree) {
                distance(obs.tree, sim.tree, config)
		    })
        }

        # SMC approximation to likelihood ratio
        old.nbhd <- sum(ws$dists[,i] < ws$epsilon)  # how many samples are in neighbourhood of data?
        new.nbhd <- sum(ws$new.dists[,i] < ws$epsilon)
        mh.ratio <- mh.ratio * new.nbhd / old.nbhd

        # accept or reject the proposal
        if (runif(1) < mh.ratio) {  # always accept if ratio > 1
            ws.accept <- ws.accept + 1
            ws$particle <- new.particle
            ws$dists[,i] <- new.dists
            ws$sim.trees[[i]] <- new.trees
        }
    }
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
    result <- list(
        niter <- 0,
        theta <- list(),
        weights <- list(),
        accept.rate <- {},
        epsilons <- {}
    )

    initialize.smc(ws)

    n.iter <- 1
    epsilon <- .Machine$double.xmax
    while (epsilon != config$final.epsilon) {
        ws$accept <- 0  # FIXME: do these have to be set again below?
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
        perturb.particles(ws)  # Metropolis-Hastings sampling happens here

        # record everything
        result$theta[[niter]] <- ws$particles
        result$weights[[ninter]] <- ws$weights
        result$epsilons <- c(result$epsilons, ws$epsilon)
        result$accept.rate <- c(result$accept.rate, ws.accept / ws.alive)

        if (!is.na(trace.file)) {
            for (i in 1:config$nparticle) {
                write.table(
                    x=c(niter, i, ws$weights[i], ws$particles[i,], ws$dists[,i]),
                    file=trace.file,
                    append=TRUE,
                    sep="\t"
                )
            }
        }

        niter <- niter + 1

        # if acceptance rate is low enough, we're done
        if (result$accept.rate[niter-1] <= config$final.accept.rate) {
            ws$epsilon <- config$final.epsilon
            break  # FIXME: this should be redundant given loop condition above
        }
    }

    # finally sample from the estimated posterior distribution
    resample.particles(ws)
    result$theta[[niter]] <- ws$particles
    result$weights[[niter]] <- ws$weights
    result$niter <- niter

    return (result)
}



