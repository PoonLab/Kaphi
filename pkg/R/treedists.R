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
# this is natively a similarity score, but turn it into a distance metric by subtracting the similarity score from the tree size
# measures the similarity between trees as the # of tips in the largest subtree identical between trees
# tree must be bifurcating
# TODO: need the tree size so you can return treesize - MAST(tree1,tree2) as the distance metric
MAST <- function(x, y) {
  # copy the trees -- this function destroys its working copy
  tree1 <- x
  tree2 <- y
  # count the tips and ensure that the trees are the same size
  numtips <- length(tree1$tip.label)
  if(numtips != length(tree2$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to comput the MAST metric') }
  
  # label nodes -- should be able to be called with `nodelabels()`
  
  # create 2D array
  numnodes <- (numtips * 2) - 1
  vals <- matrix(nrow=numnodes, ncol=numnodes)
  
  # recursively calculate mast score
  for (node1 in 1:numnodes) {
    for (node2 in 1:numnodes) {
      vals[node1,node2] <- rmast(node1, node2, tree1, tree2)
    }
  }
  
  biggest <- 0
  for (i in 1:numnodes) {
    for (j in 1:numnodes) {
      if (vals[i,j] > biggest) {
        biggest <- vals[i,j]
      }
    }
  }
  
  return(biggest)
}


# the mast score is the number of tips in the MAST of the given subtrees
rmast <- function(a, w, treea, treew) {
  # only want to narrow in on the tips present in a subtree 'rooted' with node a and a subtree 'rooted' at node w
  vector.a <- vector()
  vector.w <- vector()
  subtree.a <- .retrieve.tips(a, treea, vector.a)
  subtree.w <- .retrieve.tips(w, treew, vector.w)
  
  if (is.na(vals[a,w]) != TRUE) {   # a = node1, w = node2
    return (vals[a,w])     # already calculated
  } else {
    # leaf branch
    if (a <= length(treea$tip.label)) {
      #if (treea$tip.label[a] %in% treea$tip.label) {                # problem is comparing tip.labels, but trying to compare arbitrary integers right now
        if (treea$tip.label[a] %in% subtree.w) {              # didn't actually make an associative array...observed that tip.label vector was the same order as tiplabel index number 
          return(1)                                 # tree a
        } else {
          return(0)
        }
      #}
    }
    
    if (w <= length(treew$tip.label)) {
      #if (treew$tip.label[w] %in% treew$tip.label) {
        if (treew$tip.label[w] %in% subtree.a) {
          return(1)                                 # tree w
        } else {
          return(0)
        }
      #}
    }
    
    # non-leaf branch, this must be implemented for a bifurcating tree
    children.a <- sapply(which(a == treea$edge[,1]), function(x){treea$edge[x,2]})
    b <- min(children.a)                     # left child of node 'a' 
    c <- max(children.a)                     # right child of node 'a'
    children.w <- sapply(which(w == treew$edge[,1]), function(x){treew$edge[x,2]})
    x <- min(children.w)                     # left child of node 'w'
    y <- max(children.w)                     # right child of node 'w'
    
    step1 <- rmast(b,x,treea,treew) + rmast(c,y,treea,treew)
    step2 <- rmast(b,y,treea,treew) + rmast(c,x,treea,treew)
    step3 <- rmast(a,x,treea,treew)
    step4 <- rmast(a,y,treea,treew)
    step5 <- rmast(b,w,treea,treew)
    step6 <- rmast(c,w,treea,treew)
    results <- c(step1, step2, step3, step4, step5, step6)
    
    vals[a,w] <- max(results)
    return (vals[a,w])
  }
}

## function retrieves tip labels (NAMES) only of a subtree with a given node label
.retrieve.tips <- function(node, tree, vect) {
  children <- sapply(which(node == tree$edge[,1]), function(x){tree$edge[x,2]})
  # if a tip, record and store tip.label index
  if (length(children) == 0) {
    vect <- append(vect, tree$tip.label[node])
    return(vect)
  }
  descendants <- sapply(children, function(x){.retrieve.tips(x, tree, vect)})
  
  # clean this up
  return(unique(unlist(descendants<-append(descendants, descendants))))
}



##-------------------------------------------------------------------------------------------------------------
# Align (Nye et al., 2006)
Align <- function(x, y) {
  # make 2D array of Nye distances
  # find shortest path through them
  # copy the trees -- this function destroys its working copy
  tree1 <- x
  tree2 <- y
  
  # count the tips
  numtips <- length(tree1$tip.label)
  if (numtips != length(tree2$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to compute Align metric.') }
  t1.internals <- (numtips+1) : (numtips*2 - 1)   # why are the number of internal nodes subtracted by 2? Nnode is subtracted by 1 in R's phylo trees
  t2.internals <- t1.internals
  
  # create a 2D-matrix to store all individual scores for edge 1 to edge 2 comparison
  scores <- matrix(nrow=length(t1.internals), ncol=length(t2.internals), dimnames=list(t1.internals, t2.internals))
  
  
  ## Nye et al metric score across all combinations of edge pairs
  for (i in t1.internals) {
    for (j in t2.internals) {
      # TODO: if the internal is the "root": ie. has no ancestor
      if (i %in% tree1$edge[,2] == FALSE) {
        list.partitions <- descendant.subset(i, tree1)
        t1.left <- append(list.partitions[[1]], i)        # adding the chosen node arbitrarily to the left 
        t1.right <- list.partitions[[2]]
      } else { 
        t1.left <- append(descendant.subset(i, tree1), i)                # arbitrarily calling the descendants of the edge the left partition
        t1.right <- setdiff( c(1:(numtips*2 - 1)), t1.left)
      }  
      
      if (j %in% tree2$edge[,2] == FALSE) {
        list.partitions <- descendant.subset(j, tree2)
        t2.left <- append(list.partitions[[1]], j)
        t2.right <- list.partitions[[2]]
      } else { 
        t2.left <- append(descendant.subset(j, tree2), j)
        t2.right <- setdiff( c(1:(numtips*2 - 1)), t2.left)
      }
      
      # 4 computations, 2 on same side, 2 on opposite side
      # for each 2 sets, take the intersection / union
      val1.same.side <- length(intersect(t1.left, t2.left)) / length(union(t1.left, t2.left))   # left side of t1 w/ left side of t2
      val2.same.side <- length(intersect(t1.right, t2.right)) / length(union(t1.right, t2.right))   # right side of t1 w/ right side of t2
      val1.opp.side <- length(intersect(t1.left, t2.right)) / length(union(t1.left, t2.right))    # left side of t1 w/ right side of t2
      val2.opp.side <- length(intersect(t1.right, t2.left)) / length(union(t1.right, t2.left))    # right side of t1 w/ left side of t2
      
      # calculate score associated with edge 1 and edge 2 and store the score in a matrix
      same.side <- c(val1.same.side, val2.same.side)
      opp.side <- c(val1.opp.side, val2.opp.side)
      res <- max(c(min(same.side), min(opp.side)))
      
      # convert the Nye et al similarity to a score by subtracting from 1 (implemented in align.py from Kuhner and Yamato)
      scores[(i-numtips),(j-numtips)] <- 1.0 - res
    }
  }
  
  
  ## Munkres algorithm 
  
  
  return(scores)
}

## MUST be bifurcating trees to use this distance metric (for the moment)  Nye et al said it doesn't necessarily have to be 
## TODO: separate nodes for case of multifurcating trees
descendant.subset <- function(edge, tree) {
  edges <- tree$edge
  children <- sapply(which(edge == tree$edge[,1]), function(x){tree$edge[x,2]})
  subsetList <- vector()
  descendants  <- sapply(children, function(x) {preorder.traversal(x, subsetList, tree)})    # vector of all descendants of a given edge
  
  # for root case, return a list of left partition and right partition
  if (edge %in% tree1$edge[,2] == FALSE) {
    subset <- sapply(descendants, function(x){unique(as.vector(unlist(x, recursive=FALSE)))})
  } else {  # else, return descendants partition (ancestors partition are dealt with one level up)
    subset <- unique(unlist(descendants, recursive=FALSE))
  }
  
  return(subset)
}

## function determines descendants of a given child (node) from a tree
# for the root, separates directly into left and right partitions
# for any other edge, all the children including the given node are the resulting 'left' partition, and the 'right partition is dealt w/ one level up
preorder.traversal <- function(child, subsetList, tree) {
  if (is.null(child)) { return(subsetList) }
  
  # store the node in a list
  subsetList <- append(subsetList, child)
  
  children <- sapply(which(child == tree$edge[,1]), function(x){tree$edge[x,2]})
  if(length(children) == 0) { return(subsetList) }
  
  if (child %in% tree1$edge[,2] == FALSE) {
    # for root case, partition and store into left and right partitions
    left.child <- min(children)
    left.subset <- preorder.traversal(left.child, subsetList, tree)
    right.subset <- setdiff( c(1:(numtips*2 - 1)), left.subset)
    
    descendants <- list(c(left.subset, right.subset))
  } else {
    # traverse children of child node
    descendants <- sapply(children, function(x) {preorder.traversal(x, subsetList, tree)})
  }
  return(descendants)
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
