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

require(yaml)

load.config <- function(file) {
    # default settings
    config <- list(
        params=NA,
        priors=list(),
        prior.densities=list(),
        proposals=list(),
        proposal.densities=list(),
        model=NA,

        # SMC settings
        nparticle=1000,
        nsample=5,
        ess.tolerance=1.5,
        final.epsilon=0.01,
        final.accept.rate=0.015,
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
	for (par.name in config$params) {
		sublist <- settings$priors[[par.name]]

        rng.call <- paste('r', sublist$dist, '(n=1,', sep='')
		arguments <- sapply(sublist[['hyperparameters']], function(x) paste(names(x), x, sep='='))
		rng.call <- paste(rng.call, paste(arguments, collapse=','), ')', sep='')
		config$priors[par.name] <- rng.call

        # need to set (arg.prior) to the desired value before calling this expression
        den.call <- paste0('d', sublist$dist, '(arg.prior,', paste(arguments, collapse=','), ')')
        config$prior.densities[par.name] <- den.call
	}

    # parse proposal settings
    if (any(!is.element(names(settings$proposals), config$params))) {
        stop("Error in load.config(): variable set in proposals is not a subset of priors")
    }
    for (par.name in names(settings$proposals)) {
        sublist <- settings$proposals[[par.name]]

        rng.call <- paste0('r', sublist$dist, '(n=1,')
        arguments <- sapply(sublist[['parameters']], function(x) paste(names(x), x, sep='='))
        rng.call <- paste0(rng.call, paste(arguments, collapse=','), ')')
        config$proposals[par.name] <- rng.call

        # note this assumes that we pass an argument as 'arg.delta' in global namespace
        den.call <- paste0('d', sublist$dist, '(arg.delta,', paste(arguments, collapse=','), ')')
        config$proposal.densities[par.name] <- den.call
    }

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
    # 3. tips = if Int, the number of tips; if vector, the tip heights
	g.args <- names(formals(generator))
	if (length(g.args)<3 || any(!is.element(c('theta', 'nsim', 'tips'), g.args))) {
		stop("generator is not a Kaphi-compatible model")
	}
    # check that function has name attribute
    if (is.null(attr(generator, 'name'))) {
        stop("generator is not a Kaphi-compatible model")
    }

    config$model <- generator
    return(config)
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
        # add 'r' prefix for random deviate generation
        theta <- sapply(config$priors, function(e) eval(parse(text=e)))
        if (length(config$params) == length(theta)) {
            names(theta) <- config$params
        }
        return(theta)  # a named vector
    }
}

prior.density <- function(config, theta) {
    result <- 1.
    for (par.name in names(theta)) {
        if (is.element(par.name, names(config$prior.densities))) {
            arg.prior <- theta[par.name]
            result <- result * eval(parse(text=config$prior.densities[[par.name]]))
        }
    }
    return (result)
}


propose <- function(config, theta) {
    # draw new values from the proposal density
    delta <- sapply(names(config$proposals), function(par.name) {
        eval(parse(text=config$proposals[[par.name]]))
    })
    # apply proposal values to update parameter vector (theta)
    for (par.name in names(theta)) {
        if (is.element(par.name, names(delta))) {
            theta[par.name] <- theta[par.name] + delta[par.name]
        }
    }
    return(theta)
}


proposal.density <- function(config, theta, new.theta) {
    # sanity check
    if (length(intersect(names(theta), names(new.theta))) != length(theta)) {
        stop("Error in proposal.density(): mismatched named vectors for theta and new.theta")
    }
    result <- 1.
    for (par.name in names(theta)) {
        if (is.element(par.name, names(config$proposal.densities))) {
            arg.delta <- new.theta[par.name] - theta[par.name]
            result <- result * eval(parse(text=config$proposal.densities[[par.name]]))
        }
    }
    return (result)
}


print.smc.config <- function(config, ...) {
    cat('Kaphi SMC configuration\n\n')

    model.name <- attr(config$model, 'name')
    cat('Model:\n ', ifelse(is.null(model.name), "No model set.", model.name), '\n')

    cat('Priors:\n')
    for (i in 1:length(config$params)) {
        cat('  ', names(config$priors)[i], '\t', config$priors[[i]], '\n')
    }

    cat('Proposals:\n')
    for (par.name in config$params) {
        if (is.element(par.name, names(config$proposals))) {
            cat('  ', par.name, '\t', config$proposals[[par.name]], '\n')
        } else {
            cat('  no proposal set for ', par.name, '\n')
        }
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


plot.smc.config <- function(config, nreps=1000, numr=1, numc=1) {
    # numr = number of rows of plots to be displayed at one time
    # numc = number of columns of plots to be displayed at one time
    # display prior distributions

    # generate M x N-matrix where M is nreps and N is number of parameters
    #  containing samples from prior distributions
    y <- sapply(1:nreps, function(x) sample.priors(config))
    h <- apply(y, 1, hist)
    s <- 1        # counter
    # x11()
    par(ask=T)    # prompts user to 'Hit <Return> to see next plot'

    # FIXME: determine number of pages (indexed by i)
    for (i in 1:(numr * numc)){
      par(mfrow = c(numr, numc))     # multiple plot display option
      for (slot in 1:(numr * numc)){
        q <- quantile(y[,s], c(0.05, 0.95))  # 90% of the sample distribution from prior of s-th parameter
        plot(
          h[[s]],
          xlab=names(h)[1],
          main='Sample from prior distribution',
          xlim=q
        )  # TODO: Adjust x-axis for each plot
        s <- s + 1
      }
    }
}
