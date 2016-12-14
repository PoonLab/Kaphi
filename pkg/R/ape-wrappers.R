require(ape)

const.coalescent <- function(theta, nsim, n.tips, labels=NA, seed=NA) {
	"
	This function exists only as a prototype.
	It generates a tree under Kingman's coalescent where branch lengths are in
	units of N_e (effective population size) generations.  To compare against a
	time tree, they must be rescaled by generation time (\tau).  Thus, the key
	parameter (N_e * tau) determines the scale of the coalescent tree.
	"
    if (!is.element('Ne.tau', names(theta))) {
        stop('theta does not contain required parameter "Ne.tau"')
    }
	if (!is.na(seed)) {
		set.seed(seed)
	}
    if (all(is.na(labels))) {
        result <- lapply(1:nsim, function(x) {
            tree <- rcoal(n=n.tips, br='coalescent')
            tree$edge.length <- tree$edge.length * theta['Ne.tau']  # rescale
            tree
        })
    } else {
        result <- lapply(1:nsim, function(x) {
            tree <- rcoal(n=n.tips, tip.label=labels, br='coalescent')
            tree$edge.length <- tree$edge.length * theta['Ne.tau']  # rescale
            tree
        })
    }
	return(result)
}
attr(const.coalescent, 'name') <- "constant.coalescent"
