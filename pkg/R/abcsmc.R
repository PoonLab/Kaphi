# Based on adaptive SMC algorithm proposed by Del Moral, Pierre, Arnaud
# Doucet, and Ajay Jasra. "An adaptive sequential Monte Carlo method for
# approximate Bayesian computation." Statistics and Computing 22.5 (2012):
# 1009-1020.

require(yaml)
require(TreeSim)

# global parameters
resize.amount <- 100
bisection.max.iter <- 10000


# This S4 object stores run settings and methods for:
# (1) sampling particles from the prior distribution
# (2) simulating a tree from a user-specified model for given parameters
# All S4 object variables will be read-only.
# based on C struct smc_config in smc.h
smc.config <- setClass("smc.config", 
	slots=c(
	    params="character",  # vector of parameter names in the model
	    priors="character",	 # vector of R expressions to generate random variates
	    nparticle="numeric", # number of particles to approximate posterior
	    nsample="numeric",   # number of simulations per particle
	    ess.tolerance="numeric", # ESS below this value triggers resampling
	
	    final.epsilon="numeric", # tolerance level to end at
	    final.accept.rate="numeric", # MCMC acceptance rate stopping criterion
	    quality="numeric",       # between 0 (fast, coarse) and 1 (slow, accurate)
	    step.tolerance="numeric", # tolerance for bisection solution of next epsilon
	    
	    decay.factor="numeric",  # subset tree kernel parameters
	    rbf.variance="numeric",
	    sst.contorl="numeric",
	    sst.normalize="numeric",
	    norm.mode="character"
	),
	prototype=list(  # defaults		
		nparticle=10, # 1000
		nsample=1, # 5
		ess.tolerance=1.5,
		
		final.epsilon=0.01,
		final.accept.rate=0.015,
		quality=0.95,
		step.tolerance = 1e-5,
		
		decay.factor=0.2,
		rbf.variance=2.0,
		sst.control=1.0,
		sst.normalize=1.0,
		norm.mode="MEAN"
	)
)

# how to display contents of S4 object to user
setMethod(f='show', signature='smc.config', definition=function(object) {
	cat('Kaphi abc.smc S4 object\n\n')
	cat('Number of particles:', object@nparticle, '\n')
	cat('Number of samples per particle: ', object@nsample, '\n')
	cat('\nAnnealing parameters\n')
	cat('  ESS tolerance: '); print(as.integer(object@ess.tolerance))
	# more...
})


# This method should call a function that samples parameters from the 
# prior distributions.  It should be configured based on YAML input from
# the user.
setMethod(f='sample', signature='smc.config', 
	definition=function(obj, x, size, replace, prob) {
		if (length(obj@priors)==0) {
			cat('No prior distributions have been set; run load.priors().')
		} else {
			theta <- sapply(obj@priors, function(e) eval(parse(text=e)))
			if (length(obj@params) == length(theta)) {
				names(theta) <- obj@params
			}
			return(theta)
		}
	}
)

setGeneric(name="load.priors", def=function(object, file) {standardGeneric("load.priors")})
setMethod(f='load.priors', signature='smc.config', definition=function(object, file) {
	# YAML should be of the following format
	"
	'N':              # name of model parameter
	  'dist': 'rnorm' # name of a random generator in R{stats}
	  'hyperparameters':
	  - 'mean': 1.0
	  - 'sd':   1.0
	"
	# where PARAMETER is a string identifier for a model parameter
	#  DISTNAME is a string that corresponds to one of random generators in R{stats}
	#   e.g., "norm" corresponds to stats:rnorm()
	#  HYPERPARAM[#] is a string that 
	prior.list <- yaml.load_file(file)
	object@params <- names(prior.list)
	expressions <- {}  # call sapply(eval) on this to generate parameter vector
	for (par.name in names(prior.list)) {
		sublist <- prior.list[[par.name]]
		rng.call <- paste(sublist$dist, '(n=1,', sep='')
		arguments <- sapply(sublist[['hyperparameters']], function(x) paste(names(x), x, sep='='))
		rng.call <- paste(rng.call, paste(arguments, collapse=','), ')', sep='')
		expressions <- c(expressions, rng.call)
	}
	print(expressions)
	object@priors <- expressions
	return(object)
})
# Usage: foo <- load.priors('examples/example-priors.yaml', foo)


# This method should call a function that simulate a tree given a 
# named parameter vector.
# abc.smc class does not need to know about any details of the actual model
# other than parameter names and priors
setGeneric(name="simulate.trees", def=function(object, nsim, seed, params) {standardGeneric("simulate.trees")})
setMethod(f='simulate.trees', signature='smc.config', definition=function(object, nsim, seed, params) {cat('Simulation method has not yet been set.')})


setGeneric(name="set.model", def=function(object, model.name) {standardGeneric("set.model")})
setMethod(f='set.model', signature='smc.config', 
		 definition=function(object, generator) {
	# emulate how glm() handles 'family' argument
	if (is.character(generator)) {
		generator <- get(generator, mode="function", envir=parent.frame())
	}
	if (is.function(generator)) {
		generator <- generator()
	}
	if (is.null(generator$generator)) {
		print(generator)
		stop("'generator' not recognized")
	}
	# make sure parameter names line up
})


# Generators are wrappers for R functions that simulate trees



## the following functions operate on a smc.config instance
foobaz <- function (object) {
	print(object@priors)
}

# initialize.particles
initialize.smc <- function (object, ) {
	
}

# .next.epsilon <- function()

# .epsilon.objfun <- function()

# .ess <- function()


# .resample <- function()

# .perturb <- function()



run.smc <- function(object, config, seed, nthreads, obs.tree, trace.file) {
	# @param config: an instance of S4 object smc.config
	# @param seed: 
	
	# matrix of parameter vectors across particles
	theta <- matrix(NA, nrow=config@nparticle, ncol=config@nparam)
	
	# for proposals
	new.theta <- matrix(NA, nrow=config@nparticle, ncol=config@nparam)
	
	# store mean kernel scores?
	x <- matrix(NA, nrow=config@nsample, ncol=config@nparticle)
	new.x <- matrix(NA, nrow=config@nsample, ncol=config@nparticle)
	
	# weights?
	w <- rep(NA, times=config@nparticle)
	new.w <- rep(NA, times=config@nparticle)
	
	# space for returned values
	accept.rate <- {}
	epsilons <- {}
	
	## Step 0: sample particles from prior distribution
	
}


