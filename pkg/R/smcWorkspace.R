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

parse.input.tree <- function(obs.tree, config) {
  # validate config argument
  #required.args <- c('decay.factor', 'rbf.variance', 'sst.control', 'norm.mode')
  #missing.args <- required.args[!is.element(required.args, names(config))]

  #if (length(missing.args)>0) {
  #  warning(paste('parse.input.tree() missing argument', missing.args, "\n"))
  #    return (NULL)
  #}

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
	# ladderize, rescale, and compute self-kernel for the observed tree
  obs.tree <- .preprocess.tree(obs.tree, config)
  return (obs.tree)
}


# ape:node.height() does not have the expected function
# ape:branching.times() requires an ultrametric tree
.get.tip.heights <- function(phy) {
	n.tips <- Ntip(phy)
	tip.dists <- node.depth.edgelength(phy)[1:n.tips]
	return(max(tip.dists)-tip.dists)
}
.get.node.heights <- function(phy) {
  n.tips <- Ntip(phy)
  n.nodes <- Nnode(phy)
  tip.dists <- node.depth.edgelength(phy)[1:n.tips]
  node.dists <- node.depth.edgelength(phy)[(n.tips+1):(n.tips+n.nodes)]
  return(max(tip.dists)-node.dists)
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
    tip.heights=.get.tip.heights(obs.tree),
    tip.labels=ifelse(is.na(regex), NA, parse.labels(obs.tree$tip.label, regex)),

    # this will hold multiPhylo objects from particles
    sim.trees=lapply(1:config$nparticle, list),

    # smc.config S3 object
    config=config,

    # each particle is a vector of model parameters
    particles=matrix(NA, nrow=config$nparticle, ncol=nparams),

    # weights of particles
    weights=rep(NA, times=config$nparticle),
    new.weights=rep(NA, times=config$nparticle),

    # distances from kernel
    dists=matrix(NA, nrow=config$nsample, ncol=config$nparticle),

    epsilon=.Machine$double.xmax,  # current tolerance (could use Inf)

    accept=0,  # number of accepted proposals
    alive=0,    # number of live particles
    seed=NA
  )
  class(workspace) <- 'smc.workspace'
  return(workspace)
}

print.smc.workspace <- function(workspace, ...) {
  cat('Kaphi SMC workspace\n\n')
  cat('Target tree:\n')
  cat(summary(workspace$obs.tree))
  cat('Configuration:\n')
  cat(show(workspace$config))
}

