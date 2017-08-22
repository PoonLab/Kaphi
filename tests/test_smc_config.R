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

test.load.config <- function() {
    result <- load.config('tests/fixtures/test.yaml')
    expected <- list(
        params=c('x'),
        priors=list(x='rnorm(n=1,mean=0,sd=1)'),
        prior.densities=list(x='dnorm(arg.prior,mean=0,sd=1)'),
        constraints=NULL,
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
        dist='1*Kaphi::kernel.dist(x, y, decay.factor=0.1, rbf.variance=2.5, sst.control=1, norm.mode=NONE)',
        decay.factor=0.1,
        rbf.variance=2.5,
        sst.control=1.0,
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

test.parse.tips <- function() {
    tips <- 6
    result <- .parse.tips(tips)
    n.tips <- result$n.tips
    n.expected <- 6
    checkEquals(n.tips, n.expected)
    tip.heights <- result$tip.heights
    expected.heights <- c(0, 0, 0, 0, 0, 0)
    checkEquals(tip.heights, expected.heights)
    
    tips2 <- c(4, 3, 8, 2, 3)
    result2 <- .parse.tips(tips2)
    n.tips2 <- result2$n.tips
    n.expected2 <- 5
    checkEquals(n.tips2, n.expected2)
    tip.heights2 <- result2$tip.heights
    expected.heights2 <- c(4, 3, 8, 2, 3)
    checkEquals(tip.heights2, expected.heights2)
}

test.collect.tips <- function() {
  # testing default settings of collect.times and "after"
  t6 <- read.tree(text="(((A|2017:0.1,B|2017:0.1):0.15,C|2016:0.25):0.05,D|2015:0.3):0;")
  result <- collect.times(t6)
  expected <- c(2017,2017,2016,2015)
  checkEquals(result, expected)
  
  # testing with different delimiter, and "before"
  t7 <- read.tree(text="(((0$A:0.1,3$B:0.1):0.15,1$C:0.25):0.05,4$D:0.3):0;")
  result <- collect.times(t7, delim="$", tsample="before", fieldnum=1)
  expected <- c(-1,2,0,3)
  checkEquals(result, expected)
  
  # testing with delim = regexp
  t8 <- read.tree(text="(((A@@0:0.1,B@@5:0.1):0.15,C@@6:0.25):0.05,D@@3:0.3):0;")
  result <- collect.times(t8, delim='@.', regexp=TRUE, tsample="after", fieldnum=0)
  expected <- c(0,5,6,3)
  checkEquals(result, expected)
}