parse.input.tree <- function(obs.tree, config) {
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
	obs.tree <- rescale.tree(obs.tree, config$norm.mode)

    # cache self-kernel score
    obs.tree$kernel <- tree.kernel(obs.tree, obs.tree,
        lambda=config$decay.factor,
        sigma=config$rbf.variance,
        rho=config$sst.control,
        normalize=0,
        rescale.mode=config$norm.mode
    )
    return (obs.tree)
}


# ape:node.height() does not have the expected function
# ape:branching.times() requires an ultrametric tree
get.tip.heights <- function(phy) {
	n.tips <- Ntip(phy)
	tip.dists <- node.depth.edgelength(phy)[1:n.tips]
	return(max(tip.dists)-tip.dists)
}


parse.labels <- function(labels, regex) {
	m <- regexpr(regex, labels, perl=TRUE)
	positions <- attr(m, 'capture.start')
	lengths <- attr(m, 'capture.length')
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


init.workspace <- function(obs.tree, config, regex=NA) {
    nparams <- length(config$params)
    workspace <- list(
        # the data!
        obs.tree=parse.input.tree(obs.tree, config),  # a phylo object

        # store some useful info
        n.tips=Ntip(obs.tree),
        tip.heights=get.tip.heights(obs.tree),
        tip.labels=ifelse(is.na(regex), NA, parse.labels(obs.tree$tip.label, regex)),

        # this will hold multiPhylo objects from particles
        sim.trees=lapply(1:config$nparticle, list),

        # smc.config S3 object
        config=config,

        # each particle is a vector of model parameters
        # FIXME: should we bother to allocate these here?  see initialize.smc()
        particles=matrix(NA, nrow=config$nparticle, ncol=nparams),
        new.particles=matrix(NA, nrow=config$nparticle, ncol=nparams),

        # weights of particles
        weights=rep(NA, times=config$nparticle),
        new.weights=rep(NA, times=config$nparticle),

        # distances from kernel
        dists=matrix(NA, nrow=config$nsample, ncol=config$nparticle),
        new.dists=matrix(NA, nrow=config$nsample, ncol=config$nparticle),

        epsilon=.Machine$double.xmax,  # current tolerance (could use Inf)

        accept=0,  # number of accepted proposals
        alive=0,    # number of live particles
        seed=NA
    )
    class(workspace) <- 'smc.workspace'
    return(workspace)
}

print.smc.workspace <- function(workspace) {
    cat('Kaphi SMC workspace\n\n')
    cat('Target tree:\n', workspace$obs.tree)
}

