require(Kaphi)
require(RUnit)
source('fixtures/.R')

expected <- 1
checkEquals(expected, treekernel)
