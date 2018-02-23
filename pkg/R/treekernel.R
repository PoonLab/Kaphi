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

# DEPRECATED
# utk <- function(t1, t2, config) {
# convenience wrapper for unlabelled tree shape kernel
#  result <- tree.kernel(t1, t2,
#                        lambda=config$decay.factor,
#                        sigma=config$rbf.variance,
#                        rho=as.double(config$sst.control),
#                        normalize=0
#                        )
#  return(result)
# }


# formally 'distance'
kernel.dist <- function(t1, t2, decay.factor, rbf.variance, sst.control, rescale.mode, labelPattern, labelReplacement, gamma) {
  if (is.null(t1$kernel)) {
    stop("t1 missing self kernel in distance()")
  }
  if (is.null(t2$kernel)) {
    stop("t2 missing self kernel in distance()")
  }

    # rescale branch lengths
    nt1 <- .rescale.tree(t1, rescale.mode)
    nt2 <- .rescale.tree(t2, rescale.mode)

  k12 <- tree.kernel(
    nt1,
    nt2,
    lambda=decay.factor,
    sigma=rbf.variance,
    rho=sst.control,
    regexPattern = labelPattern,
    regexReplacement = labelReplacement,
    gamma=gamma
  )

    # we can no longer cache a tree's kernel score to itself because a distance may potentially
    # comprise more than one kernel
    k11 <- tree.kernel(
    nt1,
    nt1,
    lambda=decay.factor,
    sigma=rbf.variance,
    rho=sst.control,
    regexPattern = labelPattern,
    regexReplacement = labelReplacement,
    gamma=gamma
  )

    k22 <- tree.kernel(
    nt2,
    nt2,
    lambda=decay.factor,
    sigma=rbf.variance,
    rho=sst.control,
    regexPattern = labelPattern,
    regexReplacement = labelReplacement,
    gamma=gamma
  )

  #result <- 1. - k / sqrt(t1$kernel * t2$kernel)
    result <- 1. - k12 / sqrt(k11 * k22)
  if (result < 0 || result > 1) {
    stop(
      cat("ERROR: kernel.dist() value outside range [0,1].\n",
          "k12: ", k12, "\n",
          "k11: ", k11, "\n",
          "k22: ", k22, "\n"
      )
    )
  }
  if (is.nan(result)) {
    cat("k11:", k11, "\n")
    cat("k22:", k22, "\n")
  }
  return (result)
}


tree.kernel <- function(tree1, tree2,
                        lambda,        # decay factor
                        sigma,         # RBF variance parameter
                        rho=1.0,         # SST control parameter; 0 = subtree kernel, 1 = subset tree kernel
                        normalize=0,   # normalize kernel score by sqrt(k(t1,t1) * k(t2,t2))
                        regexPattern="",     # arguments for labeled tree kernel
                        regexReplacement="",
                        gamma=0        # label factor
                        ) {
  # make labels
  use.label <- if (any(is.na(label1)) || any(is.na(label2)) || is.null(label1) || is.null(label2)) {
    FALSE
  } else {
    new_label1 <- gsub(regexPattern, regexReplacement, tree1$tip.label)
    new_label2 <- gsub(regexPattern, regexReplacement, tree2$tip.label)
    TRUE
  }
    
  nwk1 <- .to.newick(tree1)
  nwk2 <- .to.newick(tree2)
        
#    # make labels
#    if (any(is.na(label1)) || any(is.na(label2)) || is.null(label1) || is.null(label2)) {
#        new_label1 <- new_label2 <- NA
#    } else {
#	     label <- unique(label1, label2)
#        new_label1 <- sapply(label1, function(x) which(x == label))
#        new_label2 <- sapply(label2, function(x) which(x == label))
#    }
		
  res <- .Call("R_Kaphi_kernel",
                 nwk1, nwk2, lambda, sigma, as.double(rho), use.label, gamma, normalize,
                 PACKAGE="Kaphi")
  return (res)
}
