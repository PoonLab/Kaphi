require(Kaphi)
require(RUnit)
require(yaml)

source('tests/fixtures/simple-trees.R')

test.load.config <- function() {
    result <- load.config('tests/fixtures/test.yaml')
    expected <- list(
        params=c('x'),
        priors=c('rnorm(n=1,mean=0,sd=1)'),
        model=NA,
        nparticle=10,  # all non-default values
        nsample=6,
        ess.tolerance=1.1,
        final.epsilon=0.02,
        quality=0.9,
        step.tolerance=1e-4,
        decay.factor=0.1,
        rbf.variance=2.5,
        sst.control=0.0,
        sst.normalize=0.0,
        norm.mode='NONE'
    )
    names(expected$priors) <- c('x')
    class(expected) <- 'smc.config'
    checkEquals(expected, result, checkNames=FALSE)
}
