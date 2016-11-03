
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
abc.smc <- setClass("abc.smc", 
	slots=c(
	    params="character",    # vector of parameter names in the model
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

# need to define the following class methods:
#  sample.prior() : generate parameter vector from prior distribution

setMethod(f='show', signature='abc.smc', definition=function(object) {
	cat('Kaphi abc.smc object\n')
	cat('Number of particles: '); print(as.integer(object@nparticle))
})

setMethod(f='sample', signature='abc.smc', definition=function(x, size, replace, prob) {cat('No prior distributions have been set for this abc.smc object yet.')} )

setMethod(f='simulate', signature='abc.smc', definition=function(object, nsim, seed) {cat('Simulation method has not yet been set.')})

set.priors(yaml, abc.smc.obj) {
	
}



load.priors <- function(file) {
	# YAML should be of the following format
	"
	'N':              # name of model parameter
	  'dist': 'rnorm' # name of a random generator in R{stats}
	  'mean': 1.0
	  'sd':   1.0
	"
	
	# where PARAMETER is a string identifier for a model parameter
	#  DISTNAME is a string that corresponds to one of random generators in R{stats}
	#   e.g., "norm" corresponds to stats:rnorm()
	#  HYPERPARAM[#] is a string that 
	priors <- yaml.load_file(file)
	
}


run.smc <- function(config, seed, nthreads, obs.tree, trace.file) {
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
	accept.rate <- rep(NA, times=resize.amount)
	epsilons <- rep(NA, times=resize.amount)
	
	## Step 0: sample particles from prior distribution
	
}


