require(Kaphi)
require(RUnit)
source('tests/fixtures/owls.R')  # load test fixture, g

# a production count is the degree size of internal nodes in order of vertex IDs
# the order for this tree is Node1, Node2, Node3, Strix, Asio, Athene, Tyto
# where the tree is:
#  (Tyto,(Athene,(Strix,Asio):Node3):Node2)Node1
expected <- c(2, 2, 3, 0, 0, 0, 0)
result <- get.productions(g)  # send R igraph to C and extract length attribute there

checkEquals(expected, result)
