require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.simulate.tree <- function() {
    config <- list(
        nparticle=10,
        model=const.coalescent,
        nsample=10,
        norm.mode='NONE'
    )
    workspace <- init.workspace(t2, config)
    theta <- c(1)
    names(theta) <- c('Ne.tau')
    result <- simulate.tree(workspace, theta)
    checkEquals(length(result), 10)
    checkEquals(class(result), 'list')
    st1 <- result[[1]]
    checkEquals(class(st1), 'phylo')
    checkEquals(Ntip(st1), 3)  # t2 has three tips
}
