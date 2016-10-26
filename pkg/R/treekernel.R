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
