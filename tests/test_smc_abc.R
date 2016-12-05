require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.simulate.trees <- function() {
    config <- list(
        nparticle=10,
        model=const.coalescent,
        nsample=10,
        decay.factor=0.2,
        rbf.variance=2.0,
        sst.control=1.0,
        norm.mode='NONE'
    )
    workspace <- init.workspace(t2, config)
    theta <- c(1)
    names(theta) <- c('Ne.tau')
    result <- simulate.trees(workspace, theta)
    checkEquals(length(result), 10)
    checkEquals(class(result), 'list')
    st1 <- result[[1]]

    checkEquals(class(st1), 'phylo')
    checkEquals(Ntip(st1), 3)  # t2 has three tips
    checkTrue(!is.null(st1$kernel))
}

test.ess <- function() {
    result <- ess(seq(0, 1, 0.1))
    expected <- 1./(0+0.01+0.04+0.09+0.16+0.25+0.36+0.49+0.64+0.81+1.)
    checkEquals(expected, result)
}

test.epsilon.obj.func <- function() {
    ws <- list(
        epsilon=0.17,  # previous epsilon
        dists=matrix(c(0.2,0.1,0.18,0.12,0.16,0.14), nrow=2, ncol=3),
        weights=c(0.2,0.3,0.5),
        config=list(
            nparticle=3,
            nsample=2,
            alpha=0.9
        )
    )
    result <- epsilon.obj.func(ws, 0.13)
    # particle 1: num = 1, denom = 1
    # particle 2: num = 1, denom = 1
    # particle 3: num = 0, denom = 2
    new.weights <- c(0.2 * (1/1), 0.3 * (1/1), 0.5 * 0) / 0.5
    expected <- ess(new.weights) - 0.9 * ess(ws$weights)
    checkEquals(expected, result)
}

test.next.epsilon <- function() {
    ws <- list(
        epsilon=0.11,  # previous epsilon
        dists=matrix(seq(0.01, 0.1, 0.01), nrow=1, ncol=10),
        weights=rep(0.1, times=10),
        config=list(
            nparticle=10,
            nsample=1,
            alpha=0.9,
            step.tolerance=1e-4,
            final.epsilon=0.01
        )
    )
    result <- ess(ws$weights)
    expected <- 10.
    checkEquals(expected, result)

    result <- epsilon.obj.func(ws, 0.085)  # two particles out
    expected <- 8 - 0.9 * 10
    checkEquals(expected, result)

    result <- epsilon.obj.func(ws, 0.055)  # five particles out
    expected <- 5 - 0.9 * 10
    checkEquals(expected, result)

    # adjust alpha so that root is around 0.055
    ws$epsilon <- 0.11
    ws$config$alpha <- 0.5
    result <- next.epsilon(ws)
    expected <- 0.055
    checkEquals(expected, result, tolerance=0.01)
}
