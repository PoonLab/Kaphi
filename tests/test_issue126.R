require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.issue126 <- function() {
    result <- tree.kernel(t2, t2, lambda=0.1, sigma=2.0, rho=1.0, normalize=1)
    # this is crashing with branch issue126 (vanilla igraph)
}
