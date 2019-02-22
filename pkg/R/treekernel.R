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


.rescale.tree <- function(tree, mode) {
  #print ('.rescale.tree')
  mode <- toupper(mode)
  if (!is.element(mode, c('MEAN', 'MEDIAN', 'MAX', 'NONE'))) {
    stop("Invalid mode, must be MEAN, MEDIAN, MAX or NONE")
  }
  if (mode == 'NONE') {
    return(tree)
  }
  if (mode == 'MEAN') {
    scale <- mean(tree$edge.length)
  } else if (mode == 'MEDIAN') {
    scale <- median(tree$edge.length)
  } else {
    scale <- max(tree$edge.length)
  }
  tree$edge.length <- tree$edge.length / scale
  return(tree)
}


# DEPRECATED
parse.newick <- function(tree) {
  if (class(tree)=='phylo') {
    res <- .Call("R_Kaphi_parse_newick", write.tree(tree), PACKAGE="Kaphi")
  } else if (class(tree) == 'character') {
    res <- .Call("R_Kaphi_parse_newick", tree, PACKAGE="Kaphi")
  } else {
    return (1)
  }
  return (res)
}


.preprocess.tree <- function(tree, config) {
  #print ('preprocess')
  if (class(tree) == 'character') {
    tree <- read.tree(text=tree)
  }
  if (class(tree) != 'phylo') {
    stop(".preprocess.tree() requires phylo or character (Newick) object for tree")
  }
  tree <- ladderize(tree)
  tree <- .rescale.tree(tree, config$norm.mode)
  # cache self-kernel score (only if kernel distance is desired for distance metric; specified on user-level)
  # FIXME:  this won't work for labelled kernel
  if (grepl("kernel", config$dist)) {
    tree$kernel <- tree.kernel(tree, tree,
                             lambda=config$decay.factor,
                             sigma=config$rbf.variance,
                             rho=config$sst.control,
                             normalize=0
                             )
  }
  return(tree)
}

.to.newick <- function(tree) {
  # Make sure that the tree argument is an ape phylo object
  if (class(tree)=='phylo') {
    return (write.tree(tree))
  } else if (class(tree) == 'character') {
    # make sure string is standard Newick format
    tree <- read.tree(text=tree)
    if (is.null(tree)) {
      stop(".to.newick(): String failed to parse as Newick tree string")
    }
    return (write.tree(tree))
  } else {
    stop(".to.newick(): tree argument must be a phylo or character object.")
  }
}

utk <- function(t1, t2, config) {
  # convenience wrapper for unlabelled tree shape kernel
  result <- tree.kernel(t1, t2,
                        lambda=config$decay.factor,
                        sigma=config$rbf.variance,
                        rho=as.double(config$sst.control),
                        normalize=0
                        )
  return(result)
}



# also try R-igraph?

.get.productions <- function(g) {
  # node is a tip = 0
  # internal node has two child internal nodes = 1
  # internal node has one child tip = 2
  # internal node has two child tips = 3
  n.nodes <- length(V(g))
  deg <- degree(g, mode='out')
  sapply(1:n.nodes, function(i) {
    if(deg[i]==0) {
      return(0)
    } else {
      nb <- as.integer(neighbors(g, i))
      return(sum(deg[nb] == 0)+1)
    }
  })
}

.get.branchlengths <- function(g) {
  # @return A list keyed by node name, containing either (1) a vector
  #         of branch lengths for its two children, or (2) numeric(0)
  n.nodes <- length(V(g))
  edges <- incident_edges(g, 1:n.nodes, mode='out')
  sapply(edges, function(e) get.edge.attribute(g, 'length', e))
}


.get.children <- function(g) {
  # @return A list keyed by node name, containing either:
  #         1. a vector of two child indices
  #         2. numeric(0) if this vertex is a tip (no children)
  n.nodes <- length(V(g))
  sapply(1:n.nodes, function(i) as.integer(neighbors(g, i, mode='out')))
}

.ssq <- function(x) {
  # sum of squares
  sum(x^2)
}

.tree.to.igraph <- function(tree) {
  tree <- reorder(tree, order='postorder')
  g <- as.igraph.phylo(tree)
  g <- set.edge.attribute(g, 'length', value=tree$edge.length)
  g <- set.vertex.attribute(g, 'production', value=.get.productions(g))
  g <- set.vertex.attribute(g, 'bl', value=.get.branchlengths(g))
  # store sum-of-squares for branch lengths
  g <- set.vertex.attribute(g, 'ssq.bl', 
                            value=sapply(
                              get.vertex.attribute(g, 'bl'), 
                              .ssq)
                            )
  g <- set.vertex.attribute(g, 'children', value=.get.children(g))
  g
}

.postorder <- function(g) {
  # using reorder() instead of postorder() to support `ape` version < 5.0
  # @param g: igraph object
  # @return : vector of integer indices to internal nodes in postorder traversal
  if (is.null(V(g)$production)) {
    g <- set.vertex.attribute(g, 'production', value=.get.productions(g))
  }
  which(V(g)$production > 0)
}


tree.kernel <- function(t1, t2, lambda=0.5, sigma=1.0, rho=TRUE, normalize=FALSE, 
                        label1=NA, label2=NA, gamma=0) {
  # Wrapper function for tree.kernel.R() to normalize 
  
  # convert Phylo objects to igraph
  g1 <- .tree.to.igraph(t1)
  g2 <- .tree.to.igraph(t2)
  
  k <- .tsk(g1, g2, lambda, rbf.var=sigma, sst.control=rho)
  if (normalize) {
    if (is.null(t1$kernel)) {
      t1$kernel <- .tsk(g1, g1, lambda, rbf.var=sigma, sst.control=rho)
    }
    if (is.null(t2$kernel)) {
      t2$kernel <- .tsk(g2, g2, lambda, rbf.var=sigma, sst.control=rho)
    }
    k <- k / (sqrt(t1$kernel) * sqrt(t2$kernel))
  }
  return(k)
}


.tsk <- function(g1, g2, lambda=0.5, rbf.var=1.0, sst.control=TRUE) {
  # Re-implementation of tree shape kernel in native R instead of C.
  # This will be slower but easier to maintain.
  #
  # @param g1: first tree as igraph object
  # @param g2: second tree as igraph object
  # @param lambda:  exponential decay factor to prevent large diagonal problem
  # @param rbf.var:  radial basis function variance parameter to penalize branch lengths
  # @param sst.control:  Moschitti's sigma; if FALSE, reject subset trees for subtree kernel
  
  # FIXME: this assumes that t1 and t2 are both rooted 
  # TODO: port tip label (state) handling from Python script
  
  if (rbf.var <= 0) stop("tree.kernel: rbf.var must be greater than 0")
  if (lambda <= 0 || lambda > 1) stop("tree.kernel: lambda must be within (0,1]")
  
  # store kernel scores for previous node comparisons
  cache <- matrix(NA, nrow=length(V(g1)), ncol=length(V(g2)))  
  
  k <- 0  # initialize return value
  
  # iterate over every pairing of nodes in two trees
  for (n1 in .postorder(g1)) {
    for (n2 in .postorder(g2)) {
      if (V(g1)$production[n1] == V(g2)$production[n2]) {
        # Gaussian radial basis function
        res <- lambda * exp(
          -1./rbf.var * (
            V(g1)$ssq.bl[n1] + V(g2)$ssq.bl[n2] - 
            2 * as.numeric(V(g1)$bl[[n1]] %*% V(g2)$bl[[n2]])
            )
          )
        for (cn1 in 1:2) {
          c1 <- V(g1)$children[[n1]][cn1]
          c2 <- V(g2)$children[[n2]][cn1]
          
          if (V(g1)$production[c1] != V(g2)$production[c2]) {
            res <- res * sst.control
          }
          else {
            if (V(g1)$production[c1] == 0) {
              # both terminal nodes
              res <- res * (sst.control + lambda)
            }
            else {
              # both non-terminal nodes
              if (is.na(cache[c1,c2])) {
                res <- res * sst.control
              }
              else {
                res <- res * (sst.control + cache[c1,c2])
              }
            }
          }
        }
        
        cache[n1,n2] <- res
        k <- k + res
      }
    }
  }
  
  return(k)
}

