require(Kaphi)
require(RUnit)
source('owls.R')  # load test fixture, g
checkEquals(7, node.count(g))
