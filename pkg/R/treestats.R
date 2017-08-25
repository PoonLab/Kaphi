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


##-----------------------------------------------------------------------------
## Tree Comparison Metrics from Kuhner & Yamato (2014)

# Node (Williams & Clifford, 1971)
# - just use cophenetic.phylo in ape.

# RF (Robinson & Flouds, 1981)
# - RF.dist in phangorn.

# RFL (Robinson & Flouds, 1979; Kuhner & Felsenstein, 1994)
# - KF.dist in phangorn.

# Path Distance Metric (Penny et al., 1982)
# - path.dist in phangorn.

# Trip (Critchlow et al., 1996)
Trip <- function(x, y){
  # The trees are assumed to contain identical tips (same number and labels)
  score <- 0.0
  count <- 0
  namelist <- x$tip.label
  
  # iterate through all pairs and get to MRCA
  combinations <- combn(namelist, 2)
  npairs <- length(combinations) / 2
  pairdict <- list()
  
  for (i in 1:npairs){
    pair <- combinations[,i]
    pairname <- paste(pair, collapse='')
    # most recent common ancestor (mrca) in tree x, y.
    mrcax <- getMRCA(x, pair)
    mrcay <- getMRCA(y, pair)
    # time from tips (in # edges) to mrca in tree x, y.
    #   this assumes that timefromtips in the python implementation meant
    #   the time from the node of interest to the tips level.
    timex <- node.depth(x)[mrcax] - 1
    timey <- node.depth(x)[mrcay] - 1
    # create dictionary of pairs and their mrca depths in x and y:
    #   ex. AB: 1 2
    time2mrca <- c(timex, timey)
    pairdict[[pairname]] <- time2mrca
  }
  
  # iterate through all triplets
  triples <- combn(namelist, 3)
  ntrips <- length(triples) / 3
  
  for (i in 1:ntrips){
    trip <- triples[,i]
    count <- count + 1
    pairs <- c()
    
    pairs.from.trip <- combn(trip, 2)
    npairs.from.trip <- length(pairs.from.trip) / 2
    
    for (j in 1:npairs.from.trip){
      pair <- pairs.from.trip[,j]
      pairname <- paste(pair, collapse='')
      pairs <- c(pairs, pairname)
    }
    times1 <- list()
    times2 <- list()
    for (pair in pairs){
      times1[[pair]] <- pairdict[[pair]][1]
      times2[[pair]] <- pairdict[[pair]][2]
    }
    short1 <- min(unlist(times1))
    short2 <- min(unlist(times2))
    matches <- FALSE
    for (pair in pairs){
      if (times1[[pair]] == short1 && times2[[pair]] == short2){
        matches <- TRUE
      }
    }
    if (matches == TRUE){
      score <- score + 1.0
    }
  }
  return(score/count)
}

# Triplet Length Distance (Kuhner & Yamato, 2014)
TripL <- function(x, y){
  # The trees are assumed to contain identical tips (same number and labels)
  score <- 0.0
  count <- 0
  namelist <- x$tip.label
  
  # iterate through all pairs and get to MRCA
  combinations <- combn(namelist, 2)
  npairs <- length(combinations) / 2
  pairdict <- list()
  
  for (i in 1:npairs){
    pair <- combinations[,i]
    pairname <- paste(pair, collapse='')
    # most recent common ancestor (mrca) in tree x, y.
    mrcax <- getMRCA(x, pair)
    mrcay <- getMRCA(y, pair)
    # time from tips (in # edges) to mrca in tree x, y.
    #   this assumes that timefromtips in the python implementation meant
    #   the time from the node of interest to the tips level.
    timex <- node.depth(x)[mrcax] - 1
    timey <- node.depth(x)[mrcay] - 1
    # create dictionary of pairs and their mrca depths in x and y:
    #   ex. AB: 1 2
    time2mrca <- c(timex, timey)
    pairdict[[pairname]] <- time2mrca
  }
  
  # iterate through all triplets
  triples <- combn(namelist, 3)
  ntrips <- length(triples) / 3
  
  for (i in 1:ntrips){
    trip <- triples[,i]
    count <- count + 1
    pairs <- c()
    
    pairs.from.trip <- combn(trip, 2)
    npairs.from.trip <- length(pairs.from.trip) / 2
    
    for (j in 1:npairs.from.trip){
      pair <- pairs.from.trip[,j]
      pairname <- paste(pair, collapse='')
      pairs <- c(pairs, pairname)
    }
    times1 <- list()
    times2 <- list()
    for (pair in pairs){
      times1[[pair]] <- pairdict[[pair]][1]
      times2[[pair]] <- pairdict[[pair]][2]
    }
    # x1 and x2 are the min values of times1 and times2
    x1 <- min(unlist(times1))
    x2 <- min(unlist(times2))
    
    # y1 and y2 are the time from tips for the mrca minus x1 and x2
    mrcax <- getMRCA(x, trip)
    mrcay <- getMRCA(y, trip)
    timex <- node.depth(x)[mrcax] - 1
    timey <- node.depth(x)[mrcay] - 1
    y1 <- timex - x1
    y2 <- timey - x2
    
    matches <- FALSE
    for (pair in pairs){
      if (times1[pair] == x1 && times2[pair] == x2){
        matches <- TRUE
        score <- score + ((abs(x1-x2) + abs(y1-y2)) * 2)
      }
    }
    if (matches == FALSE){
      score <- score + (x1 + x2 + y1 + y2)
    }
  }
  return(score/count)
}

# MAST (Gordon, 1980)


# Align (Nye et al., 2006)


# Sim (Hein et al., 2005)
Sim <- function(x, y){
  score = 0.0
  ab1 = 0.0
  ab2 = 0.0
  aa1 = 0.0
  aa2 = 0.0
  c1len = 0.0
  c2len = 0.0
  # find length of each tree
  
  # compute sum(a*b) for branches in both trees
  # compute sum(a*a) for banches in each tree
  
  
}
