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
    constraints=NULL,
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
    norm.mode='NONE',

    # distance settings
    # kernel, sackin, tree.width, etc
    dist=NULL,
    
    # cached kernel settings, left alone if not specified in user-provided yaml/distance string
    decay.factor=0.2,
    rbf.variance=100.0,
    sst.control=1.0
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
  
  # parse constraints (if present)
  if (!is.null(settings$constraints)) {
    elements <- strsplit(settings$constraints, ' ')  # list of param names and operators
    elements <- unlist(elements)
    con <- ''
    for (i in elements) {
      if (grepl('^[A-Za-z]+$', i)) {
        con <- paste0(con, "theta['", i, "']")
      } else {
        con <- paste0(con, i)
      }
    }
  config$constraints <- con
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
  
  # Parse & validate distance expression
  config$dist <- parse.distance(settings$distances)
  
  # Parse Kernel Settings
  #if (length(settings$distances > 1)) {
  #  if (is.element('kernel.dist', names(settings$distances))) {
  #    kernel.settings <- settings$distances[['kernel.dist']]
  #    config$decay.factor <- kernel.settings$decay.factor
  #    config$rbf.variance <- kernel.settings$rbf.variance
  #    config$sst.control <- kernel.settings$sst.control
  #  }
  #} else if (length(settings$distances == 1)) {
  #  # parse kernel settings from string
  #} 
  
  return (config)
}

parse.distance <- function(distance) {

  # Lists of accepted tree statistic functions, separated by package
  kaphi_stats <- list('kernel.dist', 'nLTT', 'sackin', 'colless', 'cophenetic', 'ladder.length', 'IL.nodes', 'tree.width', 
                      'max.delta.width', 'n.cherries', 'prop.unbalanced', 'avg.unbalance', 'pybus.gamma', 'internal.terminal.ratio')
  ape_stats <- list('balance', 'cophenetic.phylo', 'dist.nodes', 'dist.topo')
  phyloTop_stats <- list('avgLadder', 'getDepths', 'pitchforks')
  
  # Checks the method used to specify distance expression
  if (length(distance)==1 && is.character(distance)) {
    # The user has specified the distance expression as a string (format: "weight*function(args)+...")
    # List of the distance metrics
    dist.list <- strsplit(distance, "+", fixed=TRUE)[[1]]
    # Empty vector to hold parsed expressions
    dists <- c()
    
    for (d.metric in dist.list) {
      # Break the weight off of the function
      split.weight <- strsplit(d.metric, "*", fixed=TRUE)[[1]]
      weight <- split.weight[1]
      # Separate function from the arguments
      split.fn <- strsplit(split.weight[2], "(", fixed=TRUE)[[1]]
      fn <- split.fn[1]
      arguments <- split.fn[2]
      arguments <- gsub(")", "", arguments)
      # Convert argument string to list
      arguments <- strsplit(arguments, ", ", fixed=TRUE)[[1]]

      # Add package to function name
      if (is.element(fn, kaphi_stats)) {
        fn <- paste0('Kaphi::', fn)
      } else if (is.element(fn, ape_stats)) {
        fn <- paste0('ape::', fn)
      } else if (is.element(fn, phyloTop_stats)) {
        fn <- paste0('phyloTop::', fn)
      } else {
        stop(paste0(fn, ' is not a valid choice of distance metric'))
      }
      
      # Put it back together
      if (fn == 'Kaphi::kernel.dist' || fn == 'ape::dist.topo') {
        # These two functions take in two trees instead of one.
        dist.call <- paste0(weight, "*", fn, '(', paste(arguments, collapse=', '), ')')
      } else {
        # The function takes only 1 tree
        if (length(arguments) == 1) {
          # Default argument values are used
          dist.call <- paste0(weight, '*', 'abs(', fn, '(', paste(arguments, collapse=', '), 
                              ') - ', fn, '(y))')
        } else {
          # Argument values are specified
          dist.call <- paste0(weight, '*', 'abs(', fn, '(', paste(arguments, collapse=', '),
                              ') - ', fn, '(y, ', paste(arguments[2:length(arguments)], collapse=','), '))')
        }
      }
      dists <- c(dists, dist.call)
    }
    # combines vector of expressions into one string
    expression <- paste0(dists, collapse=' + ')
  
  } else {
    # The user has specified the distance expression as a YAML dictionary
    # Vector to hold each parsed distance expression
    dists <- c()
    
    for (d.metric in names(distance)){
      # sublist contains the weight and additional arguments for the function
      sublist <- distance[[d.metric]]
      
      # Add package to function name
      if (is.element(d.metric, kaphi_stats)) {
        fn <- paste0('Kaphi::', d.metric)
      } else if (is.element(d.metric, ape_stats)) {
        fn <- paste0('ape::', d.metric)
      } else if (is.element(d.metric, phyloTop_stats)) {
        fn <- paste0('ape::', d.metric)
      } else {
        stop(paste0(d.metric, ' is not a valid choice of distance metric'))
      }
      
      # Pulls the weight value from the list
      weight <- sublist$weight
      # Convert list of arguments to strings of the format: "argument=value"
      arguments <- lapply(seq_along(sublist), 
                          function(y, n, i) { paste(n[[i]], y[[i]], sep='=') }, 
                          y=sublist, n=names(sublist))
      # Drop the weight argument
      args <- arguments[2:length(arguments)]
  
      if (d.metric == 'kernel.dist' || d.metric == 'dist.topo') {
        # These two functions take in two trees instead of one.
        dist.call <- paste0(weight, "*", fn, '(x, y, ', paste(args, collapse=', '), ')')
      } else {
        # The function takes only 1 tree
        if (length(arguments) < 2) {
          dist.call <- paste0(weight, '*', 'abs(', fn, '(x) - ', fn, '(y))')
        } else {
          dist.call <- paste0(weight, '*', 'abs(', fn, '(x, ', paste(args, collapse=', '), 
                              ') - ', fn, '(y, ', paste(args, collapse=', '), '))')
        }
      }
      # Stores individual expressions
      dists <- c(dists, dist.call)
    }
    # combines vector of expressions into one string
    expression <- paste0(dists, collapse=' + ')
  }
  return(expression)
}


set.model <- function(config, generator) {
  # generator argument can either be function name or object
  if (is.character(generator)) {
    if (any(is.element(c('bisse', 'bisseness', 'bd', 'classe', 'geosse', 'musse', 'quasse', 'yule'), tolower(generator)))) {
      generator <- get('speciation.model', mode='function', envir=parent.frame())
    }
    else if (any(is.element(c('sir.nondynamic', 'sir.dynamic', 'sis', 'seir'), tolower(generator)))) {
      generator <- get('compartmental.model', mode='function', envir=parent.frame())
    }
    else if (is.element('const.coalescent', tolower(generator))) {
      generator <- get(generator, mode='function', envir=parent.frame())
    }
    else {
      stop("Not a Kaphi-compatible model.")
    }

	}
  # check that function takes the three required arguments
    # 1. theta = a named vector of parameter values (particle)
    # 2. nsim = number of simulations to generate per particle
    # 3. tips = size of tree to simulate
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
    # apply constraints to prior samples
    if (is.null(config$constraints)) {
      return(theta)
    } else {
      if (eval(parse(text=config$constraints))){  # parse condition string to expr and evaluate
        return(theta)  # a named vector
      } else {
        sample.priors(config)  # if condition is not satisfied, resample
      }
    }
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
  old.theta <- theta
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
  # apply constraints to prior samples
  if (is.null(config$constraints)) {
    return(theta)
  } else {
    if (eval(parse(text=config$constraints))){  # parse condition string to expr and evaluate
      return(theta)  # a named vector
    } else {
      propose(config, old.theta)  # if condition is not satisfied, resample
    }
  }
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
  cat('Constraints:\n')
  cat(config$constraints, '\n')
  
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

  cat('Distance settings\n')
  cat('  Expression:\t', config$dist, '\n')

}

plot.smc.config <- function(config, nreps=1000, numr=1, numc=1) {
  # numr = number of rows of plots to be displayed at one time
  # numc = number of columns of plots to be displayed at one time
  # display prior distributions
  y <- rbind(sapply(1:nreps, function(x) sample.priors(config)))
  if (nrow(y) == 1){
    rownames(y)[1] <- names(config$priors)
  }
  h <- apply(y, 1, density)
  s <- 1        # counter
  par(ask=T)    # prompts user to 'Hit <Return> to see next plot'
  n.slides <- ceiling(nrow(y) / (numr*numc)) #determine number of plot slides required according to specified dimensions
  
  for (i in 1:n.slides){
    par(mfrow = c(numr, numc))     # multiple plot display option
    for (slot in 1:(numr * numc)){
      if (s <= nrow(y)) {
        q <- quantile(y[s,], c(0.05, 0.95))   # 90% of the sample distribution from prior of s-th parameter
        plot(
          h[[s]], 
          xlab=names(h)[s], 
          main='Sample from prior distribution',
          xlim=q
        )  
        s <- s + 1
      } 
    }
  }
  par(ask=F)
}

# parse tip arguments for each model and creates either n tips of zero height if arg is an int 
# or tips of varying heights if arg is a numeric vector of non-negative values
.parse.tips <- function(tips) {
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
    tip.heights <- rep.int(0, n.tips)
  }
  tips <- list(n.tips=n.tips, tip.heights=tip.heights)
  return(tips)
}

# with a given Newick tree string representation or phylo object in R, 
# function will extract sample collection times from tip labels
# @delim user specifies a delimiter character or user can specify a regular expression 
# @fieldnum field number (corresponds to int time field)
# @tsample need to know where the collection times are situated: "before" or "after" the delimiter specified?
# @return value is a vector of node heights (where most recent tip is 0)
collect.times <- function(tree, delim="|", regexp=FALSE, tsample="after", fieldnum=0) {
  #check if class phylo / Newick object
  if (class(tree) != "phylo") stop('Tree must of class "phylo" or Newick object')
  
  if (regexp == TRUE) labels <- sapply(tree$tip.label, function(x) {strsplit(x, delim)})
  else labels <- sapply(tree$tip.label, function(x) {strsplit(x, delim, fixed=TRUE)})
  
  if (identical(tolower(tsample), "after")) {
    times <- as.integer(sapply(labels, function(x) {x[[2]]}))   # returns a vector of mode character
  } else if (identical(tolower(tsample), "before")) {
    times <- as.integer(sapply(labels, function(x) {x[[1]]}))
  } else {
    stop("Sampling times need to be before or after given delimiter/regexp")
  }
  # modify tip.heights in .parse.tips to reflect field number (number of years since X) -> nodeheights
  ptips <-.parse.tips(times)
  ndheights <- as.integer(sapply(ptips$tip.heights, function(x) {x-fieldnum}))
  
  return(ndheights)
}

