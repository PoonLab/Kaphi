require(ape)

rescale.tree <- function(tree, mode) {
    print ('rescale.tree')

    mode <- toupper(mode)
    if (!is.element(mode, c('MEAN', 'MEDIAN', 'MAX'))) {
        stop("Invalid mode, must be MEAN, MEDIAN or MAX")
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
    print ('preprocess')
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
    print ('to.newick')
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

tree.kernel <- function(tree1, tree2, lambda, sigma, rho, normalize, rescale.mode='MEAN') {
    print ('tree.kernel')
    nwk1 <- to.newick(preprocess.tree(tree1, rescale.mode))
    nwk2 <- to.newick(preprocess.tree(tree2, rescale.mode))

    res <- .Call("R_Kaphi_kernel",
                 nwk1, nwk2, lambda, sigma, rho, normalize,
                 PACKAGE="Kaphi")
    return (res)
}
