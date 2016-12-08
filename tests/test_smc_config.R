require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
require(yaml, quietly=TRUE)

source('tests/fixtures/simple-trees.R')

test.load.config <- function() {
    result <- load.config('tests/fixtures/test.yaml')
    expected <- list(
        params=c('x'),
        priors=list(x='rnorm(n=1,mean=0,sd=1)'),
        prior.densities=list(x='dnorm(arg.prior,mean=0,sd=1)'),
        proposals=list(x='rnorm(n=1,mean=0,sd=0.1)'),
        proposal.densities=list(x='dnorm(arg.delta,mean=0,sd=0.1)'),
        model=NA,
        nparticle=10,  # all non-default values
        nsample=6,
        ess.tolerance=1.1,
        final.epsilon=0.02,
        final.accept.rate=0.014,
        quality=0.9,
        step.tolerance=1e-4,
        decay.factor=0.1,
        rbf.variance=2.5,
        sst.control=0.0,
        norm.mode='NONE'
    )
    names(expected$priors) <- c('x')
    class(expected) <- 'smc.config'
    checkEquals(expected, result, checkNames=FALSE)
}


test.set.model <- function() {
    config <- load.config('tests/fixtures/test.yaml')
    checkException(set.model(config, ls))

    config <- set.model(config, const.coalescent)
    result <- attr(config$model, 'name')
    expected <- 'constant.coalescent'
    checkEquals(expected, result)
}


test.sample.priors <- function() {
    config <- load.config('tests/fixtures/test.yaml')

    theta <- sample.priors(config)
    checkEquals(names(theta), c('x'))

    samples <- sapply(1:1000, function(x) sample.priors(config))
    # standard error is sd/sqrt(n)
    checkEquals(mean(samples), 0, tolerance=0.1)
    checkEquals(sd(samples), 1, tolerance=0.1)
}

test.prior.density <- function() {
    config <- load.config('tests/fixtures/test.yaml')
    theta <- c(x=-0.5, z=-0.1)
    result <- prior.density(config, theta)
    expected <- 0.3520653  # dnorm(-0.5, mean=0, sd=1)
    checkEquals(expected, result, tol=1e-5, check.names=FALSE)
}

test.propose <- function() {
    config <- load.config('tests/fixtures/test.yaml')
    theta <- c(x=0.1, y=0.2)
    new.theta <- propose(config, theta)
    checkEquals(class(new.theta), "numeric")
    checkEquals(length(new.theta), 2)
    checkEquals(names(new.theta), c('x', 'y'))

    result <- replicate(100, propose(config, theta))

    # note RUnit relative difference = (Exp - Obs)/Exp
    checkEquals(0.1, mean(result[1,]), tolerance=0.1)  # SE = 0.1/sqrt(100)
    checkEquals(0.1, sd(result[1,]), tolerance=0.1)
    checkEquals(all(result[2,]==0.2), TRUE)  # other parameter should be untouched
}

test.proposal.density <- function() {
    config <- load.config('tests/fixtures/test.yaml')
    theta <- c(x=0.1, y=0.2)
    new.theta <- c(x=0.11, y=0.2)
    result <- proposal.density(config, theta, new.theta)
    expected <- 3.969525  # dnorm(0.01, mean=0, sd=0.1)
    checkEquals(expected, result, tol=0.001, check.names=FALSE)
}

