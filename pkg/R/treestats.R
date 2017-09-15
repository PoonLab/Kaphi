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

##-----------------------------------------------------------------------------
## Kaphi Tree Statistics

nLTT <- function(t1, t2) {
  nwk1 <- .to.newick(t1)
  nwk2 <- .to.newick(t2)
  res <- .Call("R_Kaphi_nLTT", nwk1, nwk2, PACKAGE="Kaphi")
  return(res)
}

# sum of node depths (branch lengths) from each tip to the root
sackin <- function(t1, use.branch.lengths=FALSE) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_sackin", nwk, use.branch.lengths, PACKAGE="Kaphi")
  return(res)
}

# sum of absolute differences in numbers of tips that descend from
#  left and right branches, for all internal nodes
colless <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_colless", nwk, PACKAGE="Kaphi")
  return(res)
}

cophenetic.index <- function(t1, use.branch.lengths=FALSE) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_cophenetic", nwk, use.branch.lengths, PACKAGE="Kaphi")
  return(res)
}

ladder.length <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_ladder_length", nwk, PACKAGE="Kaphi")
  return(res)
}

IL.nodes <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_il_nodes", nwk, PACKAGE="Kaphi")
  return(res)
}

tree.width <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_width", nwk, PACKAGE="Kaphi")
  return(res)
}

max.delta.width <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_max_delta_width", nwk, PACKAGE="Kaphi")
  return(res)
}

n.cherries <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_cherries", nwk, PACKAGE="Kaphi")
  return(res)
}

prop.unbalanced <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_prop_unbalanced", nwk, PACKAGE="Kaphi")
  return(res)
}

avg.unbalance <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_avg_unbalance", nwk, PACKAGE="Kaphi")
  return(res)
}

pybus.gamma <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_pybus_gamma", nwk, PACKAGE="Kaphi")
  return(res)
}

internal.terminal.ratio <- function(t1) {
  nwk <- .to.newick(t1)
  res <- .Call("R_Kaphi_internal_terminal_ratio", nwk, PACKAGE="Kaphi")
  return(res)
}


##-----------------------------------------------------------------------------
## Wrapper functions for supported metrics that output non-scalar values:

# cophenetic.phylo from ape
cophenetic.phylo.met <- function(x, y){
  matx <- cophenetic.phylo(x)
  maty <- cophenetic.phylo(y)
  if (all(rownames(matx) == rownames(maty))){
    corrcoef <- cor(c(matx), c(maty), method='kendall')
    return(corrcoef)
  } else {
    stop("cophenetic.phylo.met requires that the two trees being compared have the same tip labels")
  }
}

# dist.nodes from ape
dist.nodes.met <- function(x, y){
  matx <- cophenetic.phylo(x)
  maty <- cophenetic.phylo(y)
  if (all(rownames(matx) == rownames(maty))){
    corrcoef <- cor(c(matx), c(maty), method='kendall')
    return(corrcoef)
  } else {
    stop("dist.nodes.met requires that the two trees being compared have the same tip labels")
  }
}

# getDepths from phyloTop
getDepths.met <- function(x, type='tips'){
  # type is one of c('tips', 'nodes')
  res <- getDepths(x)
  if (type == 'tips') {
    val <- res$tipDepths
  } else if (type == 'nodes') {
    val <- res$nodeDepths
  } else {
    stop("type must be one of c('tips', 'nodes')")
  }
  return(val)
}

