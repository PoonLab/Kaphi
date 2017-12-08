##-----------------------------------------------------------------------------
## Tree Comparison Metrics from Kuhner & Yamato (2014)

# RF (Robinson & Flouds, 1981)
# - RF.dist in phangorn.

# RFL (Robinson & Flouds, 1979; Kuhner & Felsenstein, 1994)
# - KF.dist in phangorn.

# Path Distance Metric (Steel & Penny, 1993)
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
    #cat(trip, ":", matches, "\n")
  }
  return(score/count)
}



##-------------------------------------------------------------------------------------------------------------
# Triplet Length Distance (Kuhner & Yamato, 2014)
TripL <- function(x, y, k=0.65){
  # The trees are assumed to contain identical tips (same number and labels)
  if(length(x$tip.label) != length(y$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to compute the TripL metric.') }
  if(any(!is.element(x$tip.label, y$tip.label))) { stop('Tree 1 and Tree 2 must contain identical tip label sets to compute the TripL metric.')}
  score <- 0.0
  count <- 0
  namelist <- x$tip.label
  
  # iterate through all pairs and get to MRCA
  combinations <- combn(namelist, 2)
  npairs <- length(combinations) / 2
  pairdict <- list()
  
  # create a matrix and fill in all pairs
  pair.dict <- matrix(nrow=length(namelist), ncol=length(namelist), dimnames=list(namelist,namelist))
  
  for (i in 1:npairs){
    pair <- combinations[,i]
    pairname <- paste(pair, collapse='')
    # most recent common ancestor (mrca) in tree x, y.
    mrcax <- getMRCA(x, pair)
    mrcay <- getMRCA(y, pair)
    # time from tips (in # edges) to mrca in tree x, y.
    #   this assumes that timefromtips in the python implementation meant
    #   the time from the node of interest to the tips level.
    # ie. timex <- length(node depth from root to tip) SUBTRACT length(node depth from root to node)
    
    timex <- node.depth.edgelength(x)[which(x$tip.label == pair[1])] - node.depth.edgelength(x)[mrcax]
    timey <- node.depth.edgelength(y)[which(y$tip.label == pair[1])] - node.depth.edgelength(y)[mrcay]
    # create dictionary of pairs and their mrca depths in x and y:
    #   ex. AB: 1 2
    pairdict[[pairname]] <- c(timex, timey)
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
    
    
    mrcax <- getMRCA(x, trip)
    mrcay <- getMRCA(y, trip)
    # ie. timex <- length(node depth from root to tip) SUBTRACT length(node depth from root to node) 
    timex <- node.depth.edgelength(x)[which(x$tip.label == trip[1])] - node.depth.edgelength(x)[mrcax]
    timey <- node.depth.edgelength(y)[which(y$tip.label == trip[1])] - node.depth.edgelength(y)[mrcay]
    # y1 and y2 are the time from tips for the mrca minus x1 and x2
    y1 <- timex - x1
    y2 <- timey - x2
    
    matches <- FALSE
    for (pair in pairs){
      if (times1[pair] == x1 && times2[pair] == x2){
        matches <- TRUE
        score <- score + (abs(x1-x2) + abs(y1-y2)) ^ k
      }
    }
    if (matches == FALSE){
      score <- score + (x1 + x2 + y1 + y2)
    }
    #cat(trip, ":", matches, '\nMRCA x', mrcax, 'time x', timex, '\nMRCA y', mrcay, 'time y', timey, '\nscore', score, "\n\n")
  }
  #cat(count, '\n')
  return(score / count)
}





##-------------------------------------------------------------------------------------------------------------
## MAST (Gordon, 1980)
MAST <- function(tree1, tree2) {
  # count the tips and ensure that the trees are the same size
  numtips <- length(tree1$tip.label)
  if(numtips != length(tree2$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to compute the MAST metric') }
  
  # create 2D array
  numnodes <- (numtips * 2) - 1
  vals <- matrix(nrow=numnodes, ncol=numnodes)
  
  # recursively calculate mast score and store into matrix
  for (node1 in 1:numnodes) {
    for (node2 in 1:numnodes) {
      vals[node1,node2] <- .r.mast(node1, node2, tree1, tree2, vals)
    }
  }
  
  # subtract tree size from raw score associated w/ maximum identical subtrees
  raw.score <- max(vals)
  distance <- numtips - raw.score
  return(distance)
}


## calculates mast score for given subtrees from tree1 and tree2 and returns value to be stored into matrix
.r.mast <- function(a, w, treea, treew, vals) {
  # narrow down towards the tips present in subtree 'rooted' with node 'a' and subtree 'rooted' at node 'w'
  vector.a <- vector()
  vector.w <- vector()
  subtree.a <- .retrieve.tips(a, treea, vector.a)
  subtree.w <- .retrieve.tips(w, treew, vector.w)
  
  if (is.na(vals[a,w]) != TRUE) {                             # a = node1, w = node2
    return (vals[a,w])                                        # already calculated, returning value unchanged
  } else {
    # leaf branch of treea
    if (a <= length(treea$tip.label)) {
      #if (treea$tip.label[a] %in% treea$tip.label) {         # TODO: would be more robust to create an associative array between tip.label vector and tiplabel integers
        if (treea$tip.label[a] %in% subtree.w) {              # didn't actually make an associative array yet; observed that tip.label vector was the same order as tiplabel index number 
          return(1)                                 
        } else {
          return(0)
        }
      #}
    }
    # leaf branch of treew
    if (w <= length(treew$tip.label)) {
      #if (treew$tip.label[w] %in% treew$tip.label) {         # commented out b/c alternative way to get the answer as the previous line of code
        if (treew$tip.label[w] %in% subtree.a) {
          return(1)                                 
        } else {
          return(0)
        }
      #}
    }
    
    # non-leaf branch; recall that this must be implemented for a bifurcating tree
    children.a <- sapply(which(a == treea$edge[,1]), function(x){treea$edge[x,2]})
    b <- min(children.a)                     # left child of node 'a' 
    c <- max(children.a)                     # right child of node 'a'
    children.w <- sapply(which(w == treew$edge[,1]), function(x){treew$edge[x,2]})
    x <- min(children.w)                     # left child of node 'w'
    y <- max(children.w)                     # right child of node 'w'
    
    step1 <- .r.mast(b,x,treea,treew, vals) + .r.mast(c,y,treea,treew, vals)
    step2 <- .r.mast(b,y,treea,treew, vals) + .r.mast(c,x,treea,treew, vals)
    step3 <- .r.mast(a,x,treea,treew, vals)
    step4 <- .r.mast(a,y,treea,treew, vals)
    step5 <- .r.mast(b,w,treea,treew, vals)
    step6 <- .r.mast(c,w,treea,treew, vals)
    results <- c(step1, step2, step3, step4, step5, step6)
    
    val <- max(results)
    return (val)
  }
}

## function retrieves only the tip labels (NAMES) of a subtree with a given node label
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
Align <- function(tree1, tree2) {
  # make 2D array of Nye distances and find shortest path through them
  numtips <- length(tree1$tip.label)
  if (numtips != length(tree2$tip.label)) { stop('Tree 1 and Tree 2 must be the same size to compute Align metric.') }
  t1.internals <- (numtips+2) : (numtips*2 - 1)   # (numtips + 1) = root, skipped in kuhner & yamato
  t2.internals <- t1.internals
  
  # create a 2D-matrix to store all individual scores for edge 1 to edge 2 comparison
  scores <- matrix(nrow=length(t1.internals), ncol=length(t2.internals), dimnames=list(t1.internals, t2.internals))
  
  ## Nye et al metric score across all combinations of edge pairs
  for (i in t1.internals) {
    for (j in t2.internals) {
       
      t1.left <- .descendant.subset(i, tree1)                # arbitrarily calling the descendants of the edge the left partition
      t1.right <- setdiff( c(1:(numtips)), t1.left)
    
      t2.left <- .descendant.subset(j, tree2)
      t2.right <- setdiff( c(1:(numtips)), t2.left)
    
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
      scores[(i-(numtips+1)),(j-(numtips+1))] <- 1.0 - res
    }
  }
  
  ## Munkres algorithm written in pkg/R/align-munkres.R
  indM <- hungarian.alg(scores)
  distance <- sum(indM * scores)
  
  return(distance)
}

## MUST be bifurcating trees to use this distance metric (for the moment)  Nye et al said it doesn't necessarily have to be 
## TODO: separate nodes for case of multifurcating trees
.descendant.subset <- function(edge, tree) {
  edges <- tree$edge
  children <- sapply(which(edge == tree$edge[,1]), function(x){tree$edge[x,2]})
  subsetList <- vector()
  descendants  <- sapply(children, function(x) {.preorder.traversal(x, subsetList, tree)})    # vector of all descendants of a given edge
  
  subset <- unique(unlist(descendants, recursive=FALSE))
  return(subset)
}

## function determines descendants of a given child (node) from a tree
# for the root, separates directly into left and right partitions
# for any other edge, all the children including the given node are the resulting 'left' partition, and the 'right partition is dealt w/ one level up
.preorder.traversal <- function(child, subsetList, tree) {
  if (length(child) == 0) { return(subsetList) }
  # store the node in a list
  if (is.element(child, 1:length(tree$tip.label))) {
    subsetList <- c(subsetList, child)
  }
  
  children <- sapply(which(child == tree$edge[,1]), function(x){tree$edge[x,2]})
  if(length(children) == 0) { return(subsetList) }
  
  # traverse children of child node
  descendants <- sapply(children, function(x) {.preorder.traversal(x, subsetList, tree)})
  descendant.tips <- unlist(descendants)
  return(descendant.tips)
}




##-------------------------------------------------------------------------------------------------------------
# Sim (Hein et al., 2005)
Sim <- function(tree1, tree2){
  # averaging S_AB and S_BA from Hein et al. p. 213
  score <- 0.0
  ab1 <- 0.0
  ab2 <- 0.0
  aa1 <- 0.0
  aa2 <- 0.0
  
  # find total (branches) length of each tree
  c1len <- sum(tree1$edge.length)
  c2len <- sum(tree2$edge.length)
  if (c1len == 0 || c2len == 0) {
     stop("Unable to compute Sim metric because one or both of the trees has a total branch length of 0.")
  }
  
  t1 <- .find.clades.and.lengths(tree1)
  t2 <- .find.clades.and.lengths(tree2)
  t1.clades <- t1[[1]]
  t1.lengths <- t1[[2]]
  t2.clades <- t2[[1]]
  t2.lengths <- t2[[2]]
  
  # compute sum(a*b) for branches in both trees
  # compute sum(a*a) for banches in each tree
  for (i in 1:(length(tree1$tip.label)*2-1)) {
    if (i != (length(tree1$tip.label) + 1)) {
      t1.clade <- t1.clades[[i]]
      t1.length <- t1.lengths[[i]]
      aa1 <- aa1 + t1.length ^ 2
      match.ind1 <- .clade.is.present(t1.clade, t2.clades)
      if (match.ind1 > -1) {
        ab1 <- ab1 + t1.length * t2.lengths[[match.ind1]]
      }
    }
  }
  
  for (i in 1:(length(tree2$tip.label)*2-1)) {
    if (i != (length(tree2$tip.label) + 1)) {
      t2.clade <- t2.clades[[i]]
      t2.length <- t2.lengths[[i]]
      aa2 <- aa2 + t2.length ^ 2
      match.ind2 <- .clade.is.present(t2.clade, t1.clades)
      if (match.ind2 > -1) {
        ab2 <- ab2 + t2.length * t1.lengths[[match.ind2]]
      }
    }
  }
  
  aa1 <- aa1 / (c1len ^ 2)
  ab1 <- ab1 / (c1len * c2len)
  aa2 <- aa2 / (c2len ^ 2)
  ab2 <- ab2 / (c1len * c2len)

  score <- (ab1/aa1 + ab2/aa2) / 2.0
  distance <- 1.0 - score
  return(distance)
}

## helper function for Sim()
# takes a tree a constructs a list holding all possible clades and associating branch lengths for each tip and node in the tree
# @param tree = given tree to construct the clade and lengths data structure for
# @return list(clades, length)
#     * clades = clade associated with each tip and node in the tree stored in the tip/node's respective index
#     * lengths = branch length associated with each tip and node in the tree stored in the tip/node's respective index
.find.clades.and.lengths <- function(tree) {
  clades <- list()
  lengths <- list()
  
  # performs on all tips and nodes except for the root
  for (child in 1:(length(tree$tip.label)*2-1)) {
    subsetList <- vector()
    if (child <= length(tree$tip.label)) {
      clades[[child]] <- tree$tip.label[child]
    } else {
      indices <- .descendant.subset(child, tree)
      clades[[child]] <- sapply(indices, function(x) {tree$tip.label[x]})
    }
    total.height <- node.depth.edgelength(tree)[child]
    parent <- tree$edge[ which(tree$edge[,2] == child), 1]
    parent.branch.len <- node.depth.edgelength(tree)[parent]
    lengths[[child]] <- total.height - parent.branch.len
  }
  return(list(clades, lengths))
}

## helper function for Sim()
# determines if a given clade is present in a list of clades
# @param clade = single clade to be compared with given group
# @param group = list of clades to be compared with a single clade
# @return match = integer of the index in list group where matching clade is found, and returns -1 otherwise
.clade.is.present <- function(clade, group) {
  match = -1
  for (i in 1:length(group)) {
    if (setequal(clade, group[[i]])) { 
      match <- i
      break
    }
  }
  return(match)
}




##-------------------------------------------------------------------------------------------------------------
## Node distance metric derived from Williams and Clifford (1971) for k=1 AND/OR Path Distance Metric from Penny et al. (1982) for k=2
Node.dist <- function(tree1, tree2, k=0.65) {
  numtips <- length(tree1$tip.label)
  if (numtips != length(tree2$tip.label)) { stop("Tree 1 and Tree 2 must be the same size to be able to compute the Node distance") }
  t1.internals <- (numtips+2) : (numtips*2 - 1)   # (numtips + 1) = root, skipped in kuhner & yamato
  t2.internals <- t1.internals
  
  # set up 3D array for storage (leaves x leaves x tree)
  vals <- array(dim=c(2, numtips, numtips))
  
  # for each pair of tips, compare their lists and count non-shared members
  index <- 1
  for (tree in c(tree1, tree2)) {
    for (leaf1 in 1:numtips) {
      for (leaf2 in 1:numtips) {
        if (tree$tip.label[leaf1] == tree$tip.label[leaf2]) {
          vals[index, leaf1, leaf2] <- 0
        } else {
          annotList <- vector()
          leaf1.parents <- .find.ancestors(leaf1, tree, annotList)
          leaf2.parents <- .find.ancestors(leaf2, tree, annotList)
          nodepath.len <- length(setdiff(leaf1.parents, leaf2.parents)) + length(setdiff(leaf2.parents, leaf1.parents)) + 1     # common ancestor is always counted
          vals[index, leaf1, leaf2] <- nodepath.len
        }
      }
    }
    index <- index + 1
  }
  
  # subtract treeA and treeB results and sum
  total <- 0
  count <- 0
  for (i in 1:numtips) {
    for (j in 1:numtips) {
      count <- count + 1
      total <- total + abs(vals[1, i, j] - vals[2, i, j]) ^ k    # power of k
    }
  }
  
  distance <- total / as.numeric(count)
  return(distance)
}

## helper function for Node.dist()
# @param leaf = tip of the given tree
# @param tree = specified tree to find parents of given leaf
# @param annotList = empty vector to be populated with all parents of the leaf
# @return ancestors = vector of all the ancestors of the given leaf
.find.ancestors <- function(leaf, tree, annotList) {
  if (leaf == length(tree$tip.label) + 1) return(unique(annotList))
  parent <- tree$edge[ which(tree$edge[,2] == leaf), 1]
  annotList <- append(annotList, parent)
  ancestors <- .find.ancestors(parent, tree, annotList)
  return(unique(ancestors))
}
