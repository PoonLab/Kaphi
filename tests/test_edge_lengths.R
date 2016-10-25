require(Kaphi)
require(RUnit)
source('tests/fixtures/owls.R')  # load test fixture, g

expected <- E(g)$length
result <- edge.lengths(g)  # send R igraph to C and extract length attribute there

checkEquals(expected, result)
