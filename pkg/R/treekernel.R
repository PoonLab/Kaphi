require(ape)

rescale.tree <- function(tree, mode) {
    #print ('rescale.tree')
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

preprocess.tree <- function(tree, rescale.mode) {
    #print ('preprocess')
    if (class(tree) == 'character') {
        tree <- read.tree(text=tree)
    }
    if (class(tree) != 'phylo') {
        stop("preprocess.tree() requires phylo or character (Newick) object for tree")
    }
    tree.1 <- ladderize(tree)
    tree.2 <- rescale.tree(tree.1, rescale.mode)
}

to.newick <- function(tree) {
    #print ('to.newick')
    if (class(tree)=='phylo') {
        return (write.tree(tree))
    } else if (class(tree) == 'character') {
        # make sure string is standard Newick format
        tree2 <- read.tree(text=tree)
        return (write.tree(tree))
    } else {
        stop("tree argument must be a phylo or character object.")
    }
}

tree.kernel <- function(tree1, tree2, lambda, sigma, rho=1, normalize=0, label1=NA, label2=NA, gamma=0, rescale.mode='MEAN') {
    # make labels
    use.label <- if (any(is.na(label1)) || any(is.na(label2)) || is.null(label1) || is.null(label2))
        FALSE
    else {
    	tree1$tip.label <- label1
        tree2$tip.label <- label2
        TRUE
    }
    
    nwk1 <- to.newick(preprocess.tree(tree1, rescale.mode))
    nwk2 <- to.newick(preprocess.tree(tree2, rescale.mode))
        
#    # make labels
#    if (any(is.na(label1)) || any(is.na(label2)) || is.null(label1) || is.null(label2)) {
#        new_label1 <- new_label2 <- NA
#    } else {
#	     label <- unique(label1, label2)
#        new_label1 <- sapply(label1, function(x) which(x == label))
#        new_label2 <- sapply(label2, function(x) which(x == label))
#    }
		
    res <- .Call("R_Kaphi_kernel",
                 nwk1, nwk2, lambda, sigma, rho, use.label, gamma, normalize,
                 PACKAGE="Kaphi")
    return (res)
}
