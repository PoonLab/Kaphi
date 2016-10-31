require(Kaphi)
require(RUnit)
source('tests/fixtures/owls.R')  # load test fixture, g

# a production count is the degree size of internal nodes in order of vertex IDs
# the order for this tree is Node1, Node2, Node3, Strix, Asio, Athene, Tyto
# where the tree is:
#  (Tyto,(Athene,(Strix,Asio):Node3):Node2)Node1

# note these are zero-index, so Node1 is 0
expected <- c(1, 6, 2, 5, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1)
result <- get.children(g)  # send R igraph to C and extract length attribute there

checkEquals(expected, result)
