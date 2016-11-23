require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.kernel.trivial <- function() {
    expected <- 1
    result <- tree.kernel(t1, t1, lambda=0.1, sigma=2.0, rho=1, normalize=1, rescale.mode="NONE")
    # normalized kernel score of a tree against itself should always be 1
    checkEquals(expected, result)
}

test.kernel.unnormalized <- function() {
    # t1 has one internal node (n1) with two tips
    # delta(n1,n1) = lambda * RBF(sigma) * [(1 + delta(cn1, cn1)) * (1 + delta(cn2, cn2))]
    # delta(cn1, cn1) = lambda (tips)
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 1 * 1 * ((1+1) * (1+1))
    checkEquals(expected, result)

    result <- tree.kernel(t1, t1, lambda=0.5, sigma=1, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 0.5 * 1 * ((1+0.5) * (1+0.5))
    checkEquals(expected, result)

    # RBF = 1 for all sigma because branch lengths are identical
    result <- tree.kernel(t1, t1, lambda=1.0, sigma=0.5, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)

    # t1 is nested within t2
    result <- tree.kernel(t1, t2, lambda=1.0, sigma=Inf, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)

    # t3 has the same topology as t1 but different branch lengths
    result <- tree.kernel(t1, t3, lambda=1.0, sigma=Inf, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4  # no RBF penalty when sigma -> Infinity
    checkEquals(expected, result)

    # note sigma parameter for RBF *includes* the factor of two!
    result <- tree.kernel(t1, t3, lambda=1.0, sigma=1.0, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 1 * exp(-((0.3-0.1)^2 + (0.4-0.2)^2) / (1.0)) * ((1+1) * (1+1))
    checkEquals(expected, result)
}
