require(ape)
require(igraph)

test <- function(graph) {
    if (!is_igraph(graph)) {
        stop("Not an igraph object")
    }
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_test", graph, PACKAGE="Kaphi")
    res
}

s <- "owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
cat(s, file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
g <- as.igraph(tree.owls)
set.edge.attribute(g, "branch.length", value=tree.owls$edge.length)


