require(ape)

ape.coalescent <- function(object, nsim, labels, seed) {
	n.tips <- length(labels)
	set.seed(seed)
	if (all(is.na(labels))) {
		result <- rcoal(n=n.tips, br='coalescent')
	} else {
		result <- rcoal(n=n.tips, tip.label=labels, br='coalescent')
	}
	return(result)
}
