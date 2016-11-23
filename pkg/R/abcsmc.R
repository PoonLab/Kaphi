# Based on adaptive SMC algorithm proposed by Del Moral, Pierre, Arnaud
# Doucet, and Ajay Jasra. "An adaptive sequential Monte Carlo method for
# approximate Bayesian computation." Statistics and Computing 22.5 (2012):
# 1009-1020.


# global parameters
resize.amount <- 100
bisection.max.iter <- 10000


simulate.tree <- function(workspace, theta, seed=NA, ...) {
	# @param workspace: smc.workspace object
	# @param theta: parameter vector
	# @param seed: argument to set.seed()
    config <- workspace$config
    if (is.null(body(config$model))) {
        cat('Simulation method has not been set for configuration.')
    }
    result <- config$model(theta, config$nsample, workspace$n.tips, workspace$tip.heights, workspace$tip.labels, seed, ...)
    return(result)
}


initialize.smc <- function(ws, ...) {
    config <- ws$config
    nparams <- len(config$params)

    # reset workspace containers
    ws$sim.trees <- lapply(1:config$nparticle, list)
    ws$particles <- matrix(NA, nrow=config$nparticle, ncol=config@nparams)
    ws$new.particles <- matrix(NA, nrow=config$nparticle, ncol=config@nparams)
    ws$weights <- rep(NA, times=config$nparticle)
    ws$new.weights <- rep(NA, times=config$nparticle)
    ws$kscores <- matrix(NA, nrow=config$nsample, ncol=config$nparticle)
    ws$new.kscores <- matrix(NA, nrow=config$nsample, ncol=config$nparticle)

	for (i in 1:config@nparticle) {
        # sample particle from prior distribution
		ws$particles[i,] <- sample(config)

        # assign uniform weights
		ws$weights[i] <- 1./config@nparticle

		# simulate trees from particle
		ws$sim.trees[[i]] <- simulate.tree(ws, config, ...)

		# calculate kernel distances for trees
		ws$kscores[,i] <- sapply(ws$sim.trees[[i]], function(sim.tree) {
			tree.kernel(
				sim.tree,
				obs.tree,
				lambda=config$decay.factor,
				sigma=config$rbf.variance,
				rho=config$sst.control,
				normalize=config$sst.normalize,
				rescale.mode=config$norm.mode
			)
		})
	}
    cat('Initialized SMC workspace.\n')
}




# .perturb <- function()



next.epsilon <- function(config, prev.epsilon) {
	" Use bisection method to solve "

}







run.smc <- function(workspace, config, trace.file=NA, regex=NA, seed=NA, nthreads=1, ...) {
	"
	@param config: an instance of S4 object smc.config (read-only access)
	@param obs.tree: object of class 'phylo'
	@param trace.file: (optional) path to a file to write outputs
	@param seed: (optional) integer to set random seed
	@param nthreads: (optional) for running on multiple cores
	@param ...: additional arguments to pass to config@generator via simulate.tree()
	"
	# space for returned values
	accept.rate <- {}
	epsilons <- {}
	
    initialize.smc(workspace, config)

    n.iter <- 0
    epsilon <- .Machine$double.xmax
    while (epsilon != config$final.epsilon) {
        accept <- 0
        alive <- 0

        # update epsilon
        break
    }
}



