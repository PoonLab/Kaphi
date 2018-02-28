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
kernel.dist <- function(t1, t2, decay.factor, rbf.variance, sst.control, rescale.mode, labelPattern=NA, labelReplacement=NA, gamma=NA) {
  # we can no longer cache a tree's kernel score to itself because a distance may potentially
  # comprise more than one kernel
  # rescale branch lengths
  nt1 <- .rescale.tree(t1, rescale.mode)
  nt2 <- .rescale.tree(t2, rescale.mode)
  
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
  
  if (is.null(k11)) {
    stop("t1 missing self kernel in distance()")
  }
  if (is.null(k22)) {
    stop("t2 missing self kernel in distance()")
  }

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
                        lambda,                # decay factor
                        sigma,                 # RBF variance parameter
                        rho=1.0,               # SST control parameter; 0 = subtree kernel, 1 = subset tree kernel
                        normalize=0,           # normalize kernel score by sqrt(k(t1,t1) * k(t2,t2))
                        regexPattern=NA,       # how substrings that define states are extracted from tip labels
                        regexReplacement=NA,   # should be '\\1' by default to capture a single group
                        gamma=NA               # label factor, weight matrix for states of tip pairs
                        ) {
  
  # parse gamma string into matrix.. or vectors?
  splitStr <- unlist(strsplit(gamma, '}'))
  vec <- splitStr[nzchar(x=splitStr)]
  # which is states for tree1, which is states for tree2?
  if (grepl('t1', vec[1])) {
    t1.Str <- vec[1]
    t2.Str <- vec[2]
  } else {
    t1.Str <- vec[2]
    t2.Str <- vec[1]
  }
  # which has just names, which has names and values?
  if (grepl(':', t1.Str)) {
    t1.stnv <- unlist(strsplit(unlist(strsplit(t1.Str, '{', fixed=T))[2], ']'))
    t1.states <- sapply(seq_along(t1.stnv), function(x) {
      label <- unlist(strsplit(t1.stnv[x], ':', fixed=T))[1]
      gsub('[[:space:]]', '', label)
    })
    t1.labels <- t1.states[nzchar(x=t1.states)]
    
    t2.states <- unlist(strsplit(unlist(strsplit(t2.Str, '{', fixed=T))[2], ','))
    t2.labels <- sapply(seq_along(t2.states), function(x) {
      gsub('[[:space:]]', '', t2.states[x])
    })
    
    t1.values <- sapply(seq_along(t1.labels), function(x) {
      valueSet <- unlist(strsplit(t1.stnv[x], '[', fixed=T))[2]
      values <- unlist(strsplit(valueSet, ','))
      if (length(values) != length(t2.labels)) {
        stop ("Length of a single state's weight matrix values from Tree 1 does not match total number of states for Tree 2")
      }
      as.numeric(values)
    })
    
    
    gammaMat <- matrix(data=t1.values, nrow=length(t1.labels), ncol=length(t2.labels), dimnames=list(t1.labels, t2.labels))
  } else {
    t2.stnv <- unlist(strsplit(unlist(strsplit(t2.Str, '{', fixed=T))[2], ']'))
    t2.states <- sapply(seq_along(t2.stnv), function(x) {
      label <- unlist(strsplit(t2.stnv[x], ':', fixed=T))[1]
      gsub('[[:space:]]', '', label)
    })
    t2.labels <- t2.states[nzchar(x=t2.states)]
    
    t1.states <- unlist(strsplit(unlist(strsplit(t1.Str, '{', fixed=T))[2], ','))
    t1.labels <- sapply(seq_along(t1.states), function(x) {
      gsub('[[:space:]]', '', t1.states[x])
    })
    
    t2.values <- sapply(seq_along(t2.labels), function(x) {
      valueSet <- unlist(strsplit(t2.stnv[x], '[', fixed=T))[2]
      values <- unlist(strsplit(valueSet , ','))
      if (length(values) != length(t1.labels)) {
        stop ("Length of a single state's weight matrix values form Tree 2 does not match total number of states for Tree 1")
      }
      as.numeric(values)
    })
    
    gammaMat <- matrix(data=t2.values, nrow=length(t2.labels), ncol=length(t1.labels), dimnames=list(t2.labels, t1.labels))
  }
  
  # make labels
  use.label <- if (any(is.na(regexPattern)) || any(is.na(regexReplacement)) || is.null(regexPattern) || is.null(regexReplacement)) {
    FALSE
  } else {
    new_label1 <- gsub(regexPattern, regexReplacement, tree1$tip.label)
    new_label2 <- gsub(regexPattern, regexReplacement, tree2$tip.label)
    TRUE
  }
    
  nwk1 <- .to.newick(tree1)
  nwk2 <- .to.newick(tree2)
        
  res <- .Call("R_Kaphi_kernel",
                 nwk1, nwk2, lambda, sigma, as.double(rho), use.label, gamma, normalize,
                 PACKAGE="Kaphi")
  return (res)
}
