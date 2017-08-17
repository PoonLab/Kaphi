# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.

require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
require(yaml, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.simulate.trees <- function() {
    config <- list(
        nparticle=10,
        model=const.coalescent,
        nsample=10,
        dist='Kaphi::kernel.dist(x, y, decay.factor=0.2, rbf.variance=2.0, sst.control=1, norm.mode=NONE)',
        decay.factor=0.2,
        rbf.variance=2.0,
        sst.control=1.0,
        norm.mode='NONE'
    )
    workspace <- init.workspace(t2, config)
    theta <- c(1)
    names(theta) <- c('Ne.tau')
    result <- simulate.trees(workspace, theta, const.coalescent)
    checkEquals(length(result), 10)
    checkEquals(class(result), 'list')
    st1 <- result[[1]]

    checkEquals(class(st1), 'phylo')
    checkEquals(Ntip(st1), 3)  # t2 has three tips
    checkTrue(!is.null(st1$kernel))
}

test.distance <- function() {
    config <- list(
        dist='0.8*Kaphi::kernel.dist(x, y, decay.factor=0.5, rbf.variance=1.0, sst.control=1, norm.mode=NONE)',
        decay.factor=0.5,
        rbf.variance=1.0,
        sst.control=1,
        norm.mode="NONE"
    )
    checkException(result <- distance(t1, t1, config))

    t1$kernel <- tree.kernel(t1, t1,
        lambda=config$decay.factor,
        sigma=config$rbf.variance,
        rho=config$sst.control
    )
    t2$kernel <- tree.kernel(t2, t2,
        lambda=config$decay.factor,
        sigma=config$rbf.variance,
        rho=config$sst.control
    )
    result <- distance(t1, t1, config)
    expected <- 0.
    checkEquals(expected, result)

    result <- distance(t1, t2, config)
    cat(result, "\n")
}


test.ess <- function() {
    result <- .ess(seq(0, 1, 0.1))
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
            quality=0.9
        )
    )
    result <- .epsilon.obj.func(ws, 0.13)


    # particle 1: num = 1, denom = 1
    # particle 2: num = 1, denom = 1
    # particle 3: num = 0, denom = 2
    new.weights <- c(0.2 * (1/1), 0.3 * (1/1), 0.5 * 0) / 0.5
    expected <- .ess(new.weights) - 0.9 * .ess(ws$weights)
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
            quality=0.9,
            step.tolerance=1e-4,
            final.epsilon=0.01
        )
    )
    result <- .ess(ws$weights)
    expected <- 10.
    checkEquals(expected, result)

    result <- .epsilon.obj.func(ws, 0.085)  # two particles out
    expected <- 8 - 0.9 * 10
    checkEquals(expected, result)

    result <- .epsilon.obj.func(ws, 0.055)  # five particles out
    expected <- 5 - 0.9 * 10
    checkEquals(expected, result)

    # adjust alpha (quality) so that root is around 0.055
    ws$epsilon <- 0.11
    ws$config$quality <- 0.5
    result <- .next.epsilon(ws)$epsilon
    expected <- 0.055
    checkEquals(expected, result, tolerance=0.01)
}


test.initialize.smc <- function() {
    config <- load.config('tests/fixtures/test.yaml')
    ws <- init.workspace(t1, config)

    # this should fail because we haven't set the model
    checkException(initialize.smc(ws))

    config <- set.model(config, const.coalescent)
    ws <- init.workspace(t1, config)
    checkException(initialize.smc(ws))  # the model doesn't match the prior configuration

    config <- load.config('tests/fixtures/test-coalescent.yaml')
    config <- set.model(config, const.coalescent)
    ws <- init.workspace(t2, config)
    ws <- initialize.smc(ws, const.coalescent)

    # check generation of particles
    checkEquals(ws$config$nparticle, 10)
    checkEquals(nrow(ws$particles), 10)
    checkEquals(ncol(ws$particles), 1)
    checkEquals(colnames(ws$particles), c('Ne.tau'))

    # check assignment of particle weights
    checkEquals(length(ws$weights), 10)
    checkEquals(all(ws$weights==1./ws$config$nparticle), TRUE)

    # check simulation of trees from particles
    checkEquals(length(ws$sim.trees), 10)
    result <- sapply(ws$sim.trees, class)
    checkEquals(all(result=='list'), TRUE)
    result <- sapply(ws$sim.trees, length)
    checkEquals(all(result==6), TRUE)
    result <- sapply(ws$sim.trees, function(x) sapply(x, class))
    checkEquals(all(result=='phylo'), TRUE)

    # check kernel distances
    checkEquals(nrow(ws$dists), 6)  # nsample
    checkEquals(ncol(ws$dists), 10)  # nparticle
    checkEquals(all(ws$dists>=0.0), TRUE)
    checkEquals(all(ws$dists<=1.0), TRUE)
}


test.resample.particles <- function() {
    workspace <- list(
        # 2 parameters, 5 particles
        particles=matrix(seq(0.1,1.0,0.1), nrow=5, ncol=2),
        weights=seq(0.16, 0.24, 0.02),
        dists=matrix(seq(0, 0.1, length.out=25), nrow=5, ncol=5),
        config=list(nparticle=5, nsample=5)
    )
    class(workspace) <- 'smc.workspace'

    set.seed(42)  # sample(1:5,5,replace=T,prob=ws$weights) -> c(1,1,4,2,3)
    workspace <- .resample.particles(workspace)

    checkEquals(class(workspace), "smc.workspace")
    checkEquals(names(workspace), c('particles', 'weights', 'dists', 'config'))
    checkEquals(all(workspace$weights==0.2), TRUE)

    # FIXME: setting seed doesn't return expected result!
    # FIXME:   expected column indices 1,1,4,2,3
    # FIXME:   observed column indices 4,4,4,3,2
    # cat(show(workspace$particles), "\n")
    checkEquals(workspace$particles[,1], c(0.4,0.4,0.4,0.3,0.2))
}


serialize.trees <- function(trees) {
    # serialize Phylo objects in a list to facilitate comparisons
    return (paste(sapply(trees, function(x) write.tree(x))))
}


test.perturb.particles <- function() {
    config <- load.config('tests/fixtures/test-coalescent.yaml')
    config <- set.model(config, const.coalescent)

    nparticle <- config$nparticle
    theta <- c(Ne.tau=100)
    set.seed(100)
    obs.tree <- const.coalescent(theta, nsim=1, tips=20)[[1]]
    obs.tree <- parse.input.tree(obs.tree, config)

    ws <- init.workspace(obs.tree, config)
    ws <- initialize.smc(ws, const.coalescent)
    before <- ws
    set.seed(11)  # with this seed acceptance rate is 7/10
    after <- .perturb.particles(ws, const.coalescent)

    # all the particles should be alive
    checkEquals(after$alive, 10)

    # which particles changed?
    result <- sapply(1:nparticle, function(i) {
        all(before$particles[i,]==after$particles[i,])
    })
    checkEquals(sum(!result), after$accept)

    # check that the changed distances line up
    result2 <- sapply(1:nparticle, function(i) {
        all(before$dists[,i]==after$dists[,i])
    })
    checkEquals(all(result==result2), TRUE)

    # check which entries in tree list were updated
    result3 <- sapply(1:nparticle, function(i) {
        x <- serialize.trees(before$sim.trees[[i]])
        y <- serialize.trees(after$sim.trees[[i]])
        return(all(x==y))
    })
    checkEquals(all(result == result3), TRUE)
}


test.run.smc <- function() {
    DEACTIVATED("Used for debugging; produces a lot of output!")

    # just attempting the main loop to see if any exceptions get thrown
    # Ne.tau prior has mean exp(5)=148.4
    config <- load.config('tests/fixtures/test-coalescent.yaml')
    config <- set.model(config, const.coalescent)
    config$nparticle <- 20
    config$nsample <- 5

    nparticle <- config$nparticle
    theta <- c(Ne.tau=1000)
    set.seed(100)
    obs.tree <- const.coalescent(theta, nsim=1, tips=20)[[1]]

    ws <- init.workspace(obs.tree, config)
    result <- run.smc(ws, verbose=TRUE)
}
