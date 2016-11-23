require(Kaphi)
require(RUnit)
source('tests/fixtures/simple-trees.R')

test.kernel.trivial <- function() {
    expected <- 1
    result <- tree.kernel(t1, t1, lambda=0.1, sigma=2.0, rho=1, normalize=1, rescale.mode="NONE")
    # normalized kernel score of a tree against itself should always be 1
    checkEquals(expected, result)
}

test.kernel.unnormalized <- function() {
    result <- tree.kernel(t1, t1, lambda=0.1, sigma=2.0, rho=1, normalize=0, rescale.mode='NONE')
}
