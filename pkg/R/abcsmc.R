
# Based on adaptive SMC algorithm proposed by Del Moral, Pierre, Arnaud
# Doucet, and Ajay Jasra. "An adaptive sequential Monte Carlo method for
# approximate Bayesian computation." Statistics and Computing 22.5 (2012):
# 1009-1020.

require(yaml)
require(stats)


# global parameters
resize.amount <- 100
bisection.max.iter <- 10000

# based on C struct smc_config in smc.h
smc.config <- setClass("smc.config", 
	slots=c(
	    params="character",  # vector of parameter names in the model
	    priors="character",	 # vector of R expressions to generate random variates
	    generator="function",  # stores method to simulate tree
	    
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
	    norm.mode="character"   # normalize branch lengths?
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

# need to define the following class methods:
#  sample.prior() : generate parameter vector from prior distribution

setMethod(f='show', signature='smc.config', definition=function(object) {
	cat('Kaphi smc.config S4 object\n\n')
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
			cat('No prior distributions have been set for this abc.smc object yet.')
		} else {
			theta <- sapply(obj@priors, function(e) eval(parse(text=e)))
			if (length(obj@params) == length(theta)) {
				names(theta) <- obj@params
			}
			return(theta)  # a named vector
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
	object@priors <- expressions
	return(object)
})
# Usage: foo <- load.priors(foo, 'examples/example-priors.yaml')


"
This method should call a function that simulates trees given parameter vector
abc.smc class does not need to know about any details of the actual model
other than parameter names and priors.
"
setGeneric(name="simulate.tree", def=function(object, theta, n.tips, tip.heights=NA, labels=NA, seed=NA, ...) {standardGeneric("simulate.tree")})
setMethod(f='simulate.tree', signature='smc.config', definition=function(object, theta, nsim, n.tips, tip.heights=NA, labels=NA, seed=NA, ...) { 
	# @param object: smc.config S4 object
	# @param theta: parameter vector
	# @param n.tips: number of tips in each tree
	# @param tip.heights: a numeric vector of sampling times.  If not specified, then 
	#                     simulation should assume all tips sampled at same time (0).
	#                     Values greater than 0 are back in time.
	# @param labels: character vector of arbitrary order, corresponding to tips
	#                If tree is unlabeled, then use NA's.  Length = number of tips.
	# @param seed: argument to set.seed()
	if (is.null(body(object@generator))) {
		cat('Simulation method has not yet been set.')
	}
	result <- object@generator(theta, object@nsample, n.tips, labels=labels, seed=seed, ...)
	return(result)
})



"
Assign a wrapper function to S4 object variable @generator
First argument should be a named vector of parameter values.
"
setGeneric(name='set.model', def=function(object, generator) {standardGeneric('set.model')})
setMethod(f='set.model', signature='smc.config', definition=function(object, generator) {
	if (is.character(generator)) {
		generator <- get(generator, mode='function', envir=parent.frame())
	}
	# check that function takes the three required arguments
	g.args <- names(formals(generator))
	if (length(g.args)<3 || any(!is.element(c('theta', 'nsim', 'n.tips'), g.args))) {
		stop("'generator' not recognized")
	}
	object@generator <- generator
	return(object)
})


# ape:node.height() does not have the expected function
# ape:branching.times() requires an ultrametric tree
get.tip.heights <- function(phy) {
	n.tips <- Ntip(phy)
	tip.dists <- node.depth.edgelength(phy)[1:n.tips]
	return(max(tip.dists)-tip.dists)
}


parse.labels <- function(labels, regex) {
	m <- regexpr(regex, labels)
	positions <- as.vector(m)
	lengths <- attr(m, 'match.length')
	result <- {}
	for (i in 1:length(labels)) {
		if (positions[i] < 0) {
			result <- c(result, NA)
		} else {
			result <- c(result, substr(labels[i], positions[i], positions[i]+lengths[i]-1))
		}
	}
	return(result)
}


run.smc <- function(config, obs.tree, trace.file=NA, regex=NA, seed=NA, nthreads=1, ...) {
	"
	@param config: an instance of S4 object smc.config (read-only access)
	@param obs.tree: object of class 'phylo'
	@param trace.file: (optional) path to a file to write outputs
	@param seed: (optional) integer to set random seed
	@param nthreads: (optional) for running on multiple cores
	@param ...: additional arguments to pass to config@generator via simulate.tree()
	"
	# check input tree
	if (class(obs.tree)!='phylo') {
		if (class(obs.tree)=='character') {
			# attempt to parse Newick string
			obs.tree <- read.tree(text=obs.tree)
			if (is.null(obs.tree)) {
				stop("String passed for obs.tree failed to parse as Newick tree string")
			}
		}
	}
	
	# process the observed tree
	ladderize(obs.tree)
	obs.tree <- rescale.tree(obs.tree, config@norm.mode)
	n.tips <- Ntip(obs.tree)
	tip.heights <- get.tip.heights(obs.tree)
	if (is.na(regex)) {
		tip.labels <- NA
	} else {
		tip.labels <- parse.labels(obs.tree$tip.label, regex)
	}
	
	
	# the particles - a matrix of parameter vectors
	theta <- matrix(NA, nrow=config@nparticle, ncol=length(config@params))
	
	# for proposals
	new.theta <- matrix(NA, nrow=config@nparticle, ncol=length(config@params))
	
	# store kernel scores (distances) for current and proposed particles
	x <- matrix(NA, nrow=config@nsample, ncol=config@nparticle)
	new.x <- matrix(NA, nrow=config@nsample, ncol=config@nparticle)
	
	# current and proposed weights
	w <- rep(NA, times=config@nparticle)
	new.w <- rep(NA, times=config@nparticle)
	
	# space for returned values
	accept.rate <- {}
	epsilons <- {}
	
	## Step 0: sample particles from prior distribution (see smc.c:initialize())
	for (i in 1:config@nparticle) {
		theta[i,] <- sample(config)  # sample particle from prior distribution
		w[i] <- 1./config@nparticle  # assign uniform weights
		
		# simulate trees from particle
		sim.trees <- simulate.tree(config, theta[i,], n.tips, tip.heights, tip.labels, seed, ...)
		
		# calculate kernel distances for trees
		k.dists <- sapply(sim.trees, function(sim.tree) {
			tree.kernel(
				sim.tree, 
				obs.tree, 
				lambda=config@decay.factor,
				sigma=config@rbf.variance
			)
		})
	}
	
}



