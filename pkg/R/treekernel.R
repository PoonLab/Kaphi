require(ape)
require(igraph)

test <- function() {
    # call a function from test.c
    res <- .Call("R_Kaphi_test", PACKAGE="Kaphi")
    res
}

test2 <- function(x) {
    res <- .Call("R_Kaphi_test2", x, PACKAGE="Kaphi")
    res
}

node.count <- function(tree) {
    res <- .Call("R_Kaphi_nodecount", tree, PACKAGE="Kaphi")
    res
}

edge.lengths <- function(tree) {
    res <- .Call("R_Kaphi_get_edge_lengths", tree, PACKAGE="Kaphi")
    return(res)
}


rescale.tree <- function(tree, mode) {
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


get.productions <- function(tree) {
    res <- .Call("R_Kaphi_get_productions", tree, PACKAGE="Kaphi")
    return(res)
}

get.children <- function(tree) {
    res <- .Call("R_Kaphi_get_children", tree, PACKAGE="Kaphi")
    return(res)
}

get.branch.lengths <- function(tree) {
    res <- .Call("R_Kaphi_get_branch_lengths", tree, PACKAGE="Kaphi")
    return (res)
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


tree.kernel <- function(tree1, tree2, lambda, sigma, rho) {
    res <- .Call("R_Kaphi_kernel", tree1, tree2, lambda, sigma, rho, PACKAGE="Kaphi")
    return (res)
}
