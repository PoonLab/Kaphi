require(ape)

ape.coalescent <- function(theta, nsim, n.tips, labels=NA, seed=NA) {
	"
	This function exists only as a prototype - it generates a tree under 
	Kingman's coalescent where branch lengths are in units of N_e (effective 
	population size) generations.  To compare against a time tree, they must be
	rescaled by generation time (\tau).
	"
	print(farp)
	if (!is.na(seed)) {
		set.seed(seed)
	}
	if (is.na(labels)) {
		result <- rcoal(n=n.tips, br='coalescent')
	} else {
		result <- rcoal(n=n.tips, tip.label=labels, br='coalescent')
	}
	return(result)
}
