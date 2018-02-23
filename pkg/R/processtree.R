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
  #tree <- .rescale.tree(tree, config$rescale.mode)
  
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