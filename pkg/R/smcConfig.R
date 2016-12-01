require(yaml)

load.config <- function(file) {
    # default settings
    config <- list(
        params=NA,
        priors=list(),
        model=NA,

        # SMC settings
        nparticle=1000,
        nsample=5,
        ess.tolerance=1.5,
        final.epsilon=0.01,
        quality=0.95,
        step.tolerance=1e-5,

        # kernel settings
        decay.factor=0.2,
        rbf.variance=2.0,
        sst.control=1.0,
        norm.mode='MEAN'
    )
    class(config) <- 'smc.config'
    settings <- yaml.load_file(file)

    # parse prior settings
    config$params <- names(settings$priors)
    expressions <- {}
	for (par.name in config$params) {
		sublist <- settings$priors[[par.name]]
		rng.call <- paste(sublist$dist, '(n=1,', sep='')
		arguments <- sapply(sublist[['hyperparameters']], function(x) paste(names(x), x, sep='='))
		rng.call <- paste(rng.call, paste(arguments, collapse=','), ')', sep='')
		expressions <- c(expressions, rng.call)
	}
	config$priors <- expressions
    names(config$priors) <- config$params

    # FIXME: need to write a handler to interpret "1e-5" as float

    # parse SMC settings
    for (smc.set in names(settings$smc)) {
        if (!is.element(smc.set, names(config))) {
            stop("Unrecognized SMC setting in YAML: ", smc.set)
        }
        config[smc.set] <- settings$smc[smc.set]
    }

    # parse kernel settings
    for (k.set in names(settings$kernel)) {
        if (!is.element(k.set, names(config))) {
            stop("Unrecognized kernel setting in YAML: ", k.set)
        }
        config[k.set] <- settings$kernel[k.set]
    }
    return (config)
}



set.model <- function(config, generator) {
    # generator argument can either be function name or object
    if (is.character(generator)) {
		generator <- get(generator, mode='function', envir=parent.frame())
	}
	# check that function takes the three required arguments
    # 1. theta = a named vector of parameter values (particle)
    # 2. nsim = number of simulations to generate per particle
    # 3. n.tips = size of tree to simulate
	g.args <- names(formals(generator))
	if (length(g.args)<3 || any(!is.element(c('theta', 'nsim', 'n.tips'), g.args))) {
		stop("generator is not a Kaphi-compatible model")
	}
    # check that function has name attribute
    if (is.null(attr(generator, 'name'))) {
        stop("generator is not a Kaphi-compatible model")
    }

    config$model <- generator
}



#This method should call a function that samples parameters from the
#prior distributions.  It should be configured based on YAML input from
#the user.
sample.priors <- function(config) {
    if (class(config) != 'smc.config') {
        stop('ERROR: argument must be an smc.config object')
    }
    if (length(config$priors)==0) {
        cat('No prior distributions have been set.')
    } else {
        theta <- sapply(config$priors, function(e) eval(parse(text=e)))
        if (length(config$params) == length(theta)) {
            names(theta) <- config$params
        }
        return(theta)  # a named vector
    }
}



print.smc.config <- function(config) {
    cat('Kaphi SMC configuration\n\n')

    model.name <- attr(config$model, 'name')
    cat('Model:\n ', ifelse(is.null(model.name), "No model set.", model.name), '\n')

    cat('Priors:\n')
    for (i in 1:length(config$params)) {
        cat('  ', names(config$priors)[i], '\t', config$priors[i], '\n')
    }

	cat('SMC settings\n')
	cat('  Number of particles:', config$nparticle, '\n')
	cat('  Number of samples per particle:', config$nsample, '\n')
	cat('  ESS tolerance:', config$ess.tolerance, '\n')
    cat('  Final epsilon:', config$final.epsilon, '\n')
    cat('  Quality:', config$quality, '\n')
    cat('  Step tolerance:', config$step.tolerance, '\n')

    cat('Kernel settings\n')
    cat('  Decay factor:', config$decay.factor, '\n')
    cat('  RBF variance:', config$rbf.variance, '\n')
    cat('  SST control:', config$sst.control, '\n')
    cat('  Normalization:', config$norm.mode, '\n')
}


plot.smc.config <- function(config, nreps=1000) {
    # display prior distributions
    y <- sapply(1:nreps, function(x) sample.priors(config))
    h <- apply(y, 1, hist)
    plot(h[[1]], xlab=names(h)[1], main='Sample from prior distribution')  # TODO: Prompt user to go to next plot
}
