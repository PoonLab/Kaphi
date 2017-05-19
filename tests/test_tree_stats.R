require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.nLTT <- function() {
    result <- nLTT(t4, t5)
    expected <- 1/3. * 1/6.
    checkEquals(expected, result)
}

test.sackin <- function() {
    # sum of depths for all tips (node distance to root)
    result <- sackin(t1)
    expected <- 2
    checkEquals(expected, result)

    result <- sackin(t2)
    expected <- 5
    checkEquals(expected, result)

    # sum branch lengths instead of node depths
    result <- sackin(t1, TRUE)
    expected <- 0.3
    checkEquals(expected, result)

    result <- sackin(t2, TRUE)
    expected <- 0.2 + 0.3 + 0.3
    checkEquals(expected, result)
}

test.colless <- function() {
    result <- colless(t1)
    expected <- 0
    checkEquals(expected, result)

    result <- colless(t2)
    expected <- 1
    checkEquals(expected, result)
}

