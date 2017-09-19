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
    # most recent common ancestor (mrca) for the pair in tree x, y.
    mrcax <- getMRCA(x, pair)
    mrcay <- getMRCA(y, pair)
    # time from tips (in # edges) to mrca in tree x, y.
    #   this assumes that timefromtips in the python implementation meant
    #   the time from the node of interest to the tips level.
    timex <- node.depth(x)[mrcax] - 1
    timey <- node.depth(y)[mrcay] - 1
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
    if (matches == FALSE){
      score <- score + 1.0
    }
    cat(trip, ":", matches, "\n")
  }
  return(score/count)
}



##-------------------------------------------------------------------------------------------------------------

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
    timey <- node.depth(y)[mrcay] - 1
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
    timey <- node.depth(y)[mrcay] - 1
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
    cat(trip, ":", matches, "\n")
  }
  return(score/count)
}



##-------------------------------------------------------------------------------------------------------------
# MAST (Gordon, 1980)




##-------------------------------------------------------------------------------------------------------------
# Align (Nye et al., 2006)
Align <- function(x, y) {
  # make 2D array of Nye distances
  # find shortest path through them
  # copy the trees -- this function destroys its working copy
  tree1 <- x
  tree2 <- y
  
  # count the tips
  numtips <- len(tree1$tip.label)
  if (numtips != len(tree2$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to compute Align metric.') }
  internals <- numtips - 2    # why are the number of interanl nodes subtracted by 2? Nnode is subtracted by 1 in R's phylo trees
  
  
  ## label nodes
  # first, relabel tips as 0 to (numtips-1). This has to be done all at once as there might already have been tips with those names
  tipnames <- tree1$tip.label
  tiplist <- list()
  tipno <- 0
  for (tip in tipnames) {
    tiplist[tip] <- tipno
    tipno <- tipno + 1
  }
  for (leaf in tree1$tip.label){
    tree1$tip.label[leaf] <- tiplist[[leaf]]
  }
  for (leaf in tree2$tip.label){
    tree2$tip.label[leaf] <- tiplist[[leaf]]
  }
  
  # now label interal nodes, tarting with tips+1
  nodeno <- numtips + 1          # for some reason in their code they did NOT add 1 to numtips
  edges <- tree1$edge
  for (node in edges){
    append(tree1$node.label, nodeno)
  }
  
}



## post order tree traversal
postorder.traversal <- function(node, root, tree) {
  indices <- which(node == tree$edge[,1])                   # track down the node in column 1
  options <- sapply(indices, function(x){tree$edge[x,2]})   # record the descendant nodes of given node/edge
  
  # if one of the tips' descendants == NULL, return
  
  if (node == root) {
    node <- root
    
    # traverse down the 'left' side, the side that has the lowest node number in the recorded descendant nodes
    
    postorder.traversal(left.node, node, tree$edge)
    
    # traverse down the 'right' side, with the highest node number in the recorded descendant nodes
    right.node <- max(options)
    postorder.traversal(right.node, node, tree$edge)
    
    # label the root
    append(tree$node.label, node)
  }
}





##-------------------------------------------------------------------------------------------------------------
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
