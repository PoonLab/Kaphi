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

    # distance settings
    # kernel, sackin, tree.width, etc
    dist="kernel.dist(x, y, decay.factor=0.2, rbf.variance=100.0, sst.control=1.0)",
    
    # cached kernel settings, left alone if not specified in user-provided yaml/distance string
    decay.factor=0.2,
    rbf.variance=100.0,
    sst.control=1.0,
    norm.mode='NONE',
    labelPattern="",
    labelReplacement="",
    gamma=0
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
  if (settings$constraints != '') {
    elements <- strsplit(settings$constraints, ' ')  # list of param names and operators
    elements <- unlist(elements)
    con <- ''
    for (i in elements) {
      if (grepl('^[A-Za-z]+', i)) {
        # changes param to theta['param']
        con <- paste0(con, "theta['", i, "']")
      } else {
        # inserts an operator
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
  if (is.list(settings$distances)) {
    if (is.element('kernel.dist', names(settings$distances))) {
      kernel.settings <- settings$distances[['kernel.dist']]
      config$decay.factor <- kernel.settings$decay.factor
      config$rbf.variance <- kernel.settings$rbf.variance
      config$sst.control <- kernel.settings$sst.control
      config$norm.mode <- kernel.settings$norm.mode
      config$labelPattern <- kernel.settings$labelPattern
      config$labelReplacement <- kernel.settings$labelReplacement
      config$gamma <- kernel.settings$gamma
    }
  } else if (is.character(settings$distances)) {
    # parse kernel settings from string
    dist.list <- strsplit(settings$distances, "+", fixed=TRUE)[[1]]
    for (dist in dist.list) {
      if (grepl("kernel.dist", dist)) {
        match <- regexpr("\\(.+\\)", dist, perl=TRUE)
        args <- regmatches(dist, match)
        args <- gsub("[( )]", "", args)
        kernel.settings <- strsplit(args, ",", fixed=TRUE)[[1]]
        names <- c()
        values <- c()
        for (parm in kernel.settings) {
          split <- strsplit(parm, "=", fixed=TRUE)[[1]]
          name <- split[1]
          value <- split[2]
          names <- c(names, name)
          values <- c(values, value)
        }
        names(values) <- names
        config$decay.factor <- as.numeric(values["decay.factor"])
        config$rbf.variance <- as.numeric(values["rbf.variance"])
        config$sst.control <- as.numeric(values["sst.control"])
        config$norm.mode <- values["norm.mode"]
        config$labelPattern <- values["labelPattern"]
        config$labelReplacement <- values["labelReplacement"]
        config$gamma <- as.numeric(values["gamma"])
      }
    }
  } 
  return (config)
}

parse.distance <- function(distance) {
  # generate matrix of accepted tree statistic functions from 'metrics' list which can be added to over time without altering the rest of the function
  # if value of 1, only one variable required in distance function call (ie. sackin(x) - sackin(y))
  # if value of 2, two variables required in distance function call (ie. kernel.dist(x, y))
  metrics <- list( c('Kaphi::','kernel.dist',2),      
                   c('Kaphi::','nLTT',1), 
                   c('Kaphi::','sackin',1), 
                   c('Kaphi::','colless',1), 
                   c('Kaphi::','cophenetic.index',1), 
                   c('Kaphi::','ladder.length',1), 
                   c('Kaphi::','IL.nodes',1), 
                   c('Kaphi::','tree.width',1), 
                   c('Kaphi::','max.delta.width',1), 
                   c('Kaphi::','n.cherries',1), 
                   c('Kaphi::','prop.unbalanced',1), 
                   c('Kaphi::','avg.unbalance',1), 
                   c('Kaphi::','pybus.gamma',1), 
                   c('Kaphi::','internal.terminal.ratio',1), 
                   c('Kaphi::','cophenetic.phylo.met',1), 
                   c('Kaphi::','dist.nodes.met',1), 
                   c('Kaphi::','getDepths.met',1),
                   c('Kaphi::','Trip',2),
                   c('Kaphi::','TripL',2),
                   c('Kaphi::','MAST',2),
                   c('Kaphi::','Align',2),
                   c('Kaphi::','Sim',2),
                   c('Kaphi::','Node.dist',2),
                   
                   c('ape::','dist.topo',2),
                   
                   c('phyloTop::','avgLadder',1),
                   c('phyloTop::','pitchforks',1),
                   
                   c('phangorn::','RF.dist',2),
                   c('phangorn::','KF.dist',2),
                   c('phangorn::','path.dist',2)
                   
                   )
  mat <- matrix(nrow=length(metrics), ncol=3, dimnames=list(NULL,c('pkg', 'metric', 'no.vars')))
  stats <- t(sapply(seq_along(metrics), function(x) {mat[x,] <- metrics[[x]]}))                   # matrix of tree stats
  
  # Checks the method used to specify distance expression
  if (is.character(distance)) {                                                                # The user has specified the distance expression as a string
    dist.list <- strsplit(distance, '+', fixed=TRUE)[[1]]                                      # List of the distance metrics
    dists <- c()                                                                               # Empty vector to hold parsed expressions
    
    for (expr in dist.list) {
      pop.weight <- strsplit(expr, '*', fixed=TRUE)[[1]]                                       # Break the weight off of the function
      weight <- pop.weight[1]
      
      pop.function <- strsplit(pop.weight[2], '(', fixed=TRUE)[[1]]                            # Separate function name from the arguments 
      d.metric <- pop.function[1]
      arguments <- pop.function[2]
      arguments <- gsub(')', '', arguments)
      
      ind <- which(stats[,2] == d.metric)
      indiv.expr <- .generate.dist.expr(weight, d.metric, arguments, stats[ind,])
      dists <- append(dists, indiv.expr)
    }
  } else {
    dists <- c()  # Vector to hold each parsed distance expression                              # The user has specified the distance expression as a YAML dictionary
                                                                                                
    for (d.metric in names(distance)){ 
      sublist <- distance[[d.metric]]                                                           # sublist contains the weight and additional arguments for the function
      weight <- sublist$weight
      
      # strip out weight, package arguments, and convert list of arguments to strings of the format: "argument=value"
      all.args <- lapply(seq_along(sublist), function(y, n, i) { paste(n[[i]], y[[i]], sep='=') }, y=sublist, n=names(sublist))
      strip.args <- sapply(all.args, function(x) {which(!grepl('weight', x) && !grepl('package', x))})
      arguments <- all.args[ which(strip.args == 1) ]
      
      ind <- which(stats[,2] == d.metric)
      indiv.expr <- .generate.dist.expr(weight, d.metric, arguments, stats[ind,])
      dists <- append(dists, indiv.expr)                                                                      # Stores individual expressions
    }
  }
  expression <- paste0(dists, collapse=' + ')                                                     # combines vector of expressions into one string
  return(expression)
}


.generate.dist.expr <- function(weight, d.metric, arguments, statistic) {
  # Check if distance metric exists in stats matrix, then add package name to function name
  if (is.element(d.metric, statistic)) {
    fn <- paste0(statistic[1], d.metric)
  } else {
    stop(paste0(d.metric, ' is not a valid choice of distance metric'))
  }
  
  # Put expression together
  if (as.numeric(statistic[3]) == 2) {
    if (length(arguments) != 0 && !is.na(arguments)){                                                                 # metric takes in two variables in its function call instead of one.
      dist.call <- paste0(weight, '*', fn, '(x,y,', paste(arguments, collapse=','), ')')
    } else {                                                                                # If only default parameters are being used
      dist.call <- paste0(weight, '*', fn, "(x,y)")
    }
  } else {
    if (length(arguments) != 0 && !is.na(arguments)) {                                                                 # metric takes in one variable in its function call
      dist.call <- paste0(weight, '*', 'abs(', fn, '(x,', arguments, 
                          ')-', fn, '(y,', arguments, '))')
    } else { 
      # If only default parameters are being used
      dist.call <- paste0(weight, '*', 'abs(', fn, '(x)-', fn, '(y))')
    }
  }
  dist.call
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
    else if (is.element('epidemic', tolower(generator))) {
      generator <- get('epidem.model', mode='function', envir=parent.frame())
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

  cat('Kernel settings\n')
  cat('  Decay Factor:', config$decay.factor, '\n')
  cat('  RBF Variance:', config$rbf.variance, '\n')
  cat('  SST Control:', config$sst.control, '\n')
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

