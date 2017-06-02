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
require(diversitree, quietly=TRUE)

## Not sure if required:
## source('tests/fixtures/simple-trees.R')

test.speciation.model <- function() {
    ## Test BiSSE Model
    config <- load.config('tests/fixtures/test-bisse.yaml')
    config <- set.model(config, 'bisse')
    theta <- sample.priors(config)
    result <- speciation.model(theta nsim=1, tips=20, 'bisse')
    expected <- tree.bisse(c(0.1, 0.2, 0.003, 0.003, 0.01, 0.01), max.taxa=20)
    checkEquals(expected, result)

   ## Test BiSSE-ness Model
   config <- load.config('tests/fixtures/test-bisseness.yaml')
   config <- set.model(config, 'bisseness')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'bisseness')
   expected <- tree.bisseness(c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01, 0.05, 0.1,
                                0.05, 0.1), max.taxa=20)
   checkEquals(expected, result)

   ## Test Birth-Death Model
   config <- load.config('tests/fixtures/test-bd.yaml')
   config <- set.model(config, 'bd')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'bd')
   expected <- tree.bd(c(0.1, 0.003), max.taxa=20)
   checkEquals(expected, result)

   ## Test ClaSSE Model
   config <- load.config('tests/fixtures/test-classe.yaml')
   config <- set.model(config, 'classe')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'classe')
   expected <- tree.classe(c(2.5, 0.5, 0, 0, 0, 5, 2.41, 5.24, 0.5, 0),
                             max.taxa=20)
   checkEquals(expected, result)

   ## Test GeoSSE Model
   config <- load.config('tests/fixtures/test-geosse.yaml')
   config <- set.model(config, 'geosse')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'geosse')
   expected <- tree.geosse(c(1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), max.taxa=20)
   checkEquals(expected, result)

   ## Test MuSSE Model
   config <- load.config('tests/fixtures/test-musse.yaml') #TODO: fill out musse YAML
   config <- set.model(config, 'musse')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'musse')
   expected <- tree.musse(c(...), max.taxa=20)
   checkEquals(expected, result)

   ## Test QuaSSE Model TODO: Figure out how to write yaml file for this case
   config <- load.config('tests/fixtures/test-quasse.yaml')
   config <- set.model(config, 'quasse')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'quasse')
   lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
   mu <- function(x) constant.x(x, 0.03)
   char <- make.brownian.with.drift(0, 0.025)
   expected <- tree.quasse(c(lambda, mu, char), max.taxa=20)
   checkEquals(expected, result)

   ## Test Yule Models
   config <- load.config('tests/fixtures/test-yule.yaml')
   config <- set.model(config, 'yule')
   theta <- sample.priors(config)
   result <- speciation.model(theta nsim=1, tips=20, 'yule')
   expected <- tree.yule(c(0.1), max.taxa=20)
   checkEquals(expected, result)

   #TODO: test for varying nsim, tips, error messages(9),
}
