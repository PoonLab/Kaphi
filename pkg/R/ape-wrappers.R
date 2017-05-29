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

require(ape)

const.coalescent <- function(theta, nsim, tips, labels=NA, seed=NA) {
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
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
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
