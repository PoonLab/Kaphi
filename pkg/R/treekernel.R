require(ape)
require(igraph)

test <- function() {
    # call a function from test.c
    res <- .Call("R_Kaphi_test", PACKAGE="Kaphi")
    res
}

