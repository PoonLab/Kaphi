require(Kaphi)

s <- "owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
tree.owls <- read.tree(text=s)
g <- as.igraph(tree.owls)
set.edge.attribute(g, "branch.length", value=tree.owls$edge.length)

test()
