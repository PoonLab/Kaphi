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
require(diversitree, quietly=TRUE)
require(igraph)

## TODO: Call to igraph C library function (igraph_shortest_paths)
##       instead of adding igraph R package as a dependency.

source('tests/fixtures/simple-trees.R')

test.yule.model <- function() {
  ## Simulate tree under Yule model, n=20
  config <- load.config('tests/fixtures/test-yule.yaml')
  config <- set.model(config, 'yule')
  theta <- sample.priors(config)
  tree <- speciation.model(theta, nsim=1, tips=20, 'yule')[[1]]
  
  ## Use igraph function to generate matrix of pathlengths between nodes
  igtree <- as.igraph(tree)
  paths <- shortest.paths(igtree)
  
  ## Must remove paths involving internal nodes
  leaf.path <- paths[c(!grepl('nd', colnames(paths))), 
                     c(!grepl('nd', rownames(paths)))]

  ## Probability Distribution for probability of each path
  pd <- as.data.frame(table(leaf.path))
  pd <- pd[!(grepl('^0$', pd$leaf.path, perl=TRUE)), ]
  length <- length(pd$Freq)
  total <- sum(pd$Freq)
  
  ## Add column of probabilities
  prob <- c()
  for (i in seq(1:length)) {
    prob <- c(prob, as.numeric(pd$Freq[i]) / total)
  }
  pd$prob <- prob
  
  ## Add column of weighted path lengths
  w.mean <- c()
  for (i in seq(1:length)) {
    w.mean <- c(w.mean, as.numeric(pd$leaf.path[i]) * as.numeric(pd$prob[i]))
  }
  pd$w.mean <- w.mean
  
  exp <- sum(w.mean) # expected value
  # compare to analytical solution for distance between 2 leaves (4.2 in Yule paper) for n=20.
  # replace for loops with functionals for better run time
}


test.speciation.model <- function() {
    ## Test BiSSE Model
    config <- load.config('tests/fixtures/test-bisse.yaml')
    config <- set.model(config, 'bisse')
    theta <- sample.priors(config)
    
    trees1 <- sapply(50, function(x)) {
      
    tree.bisse(c(0.1, 0.2, 0.003, 0.003, 0.01, 0.01), max.taxa=20)
      
    trees2 <- speciation.model(theta, nsim=50, tips=20, 'bisse') 
      
    dists <- sapply(sim.trees, function(st) {
      pt <- parse.input.tree(st, config)
      distance(p.obs.tree, pt, config)
    })
    result <- mean(dists) 
    checkTrue(result < 0.2) # 0.2 is the max mean distance accepted

   ## Test BiSSE-ness Model
   #config <- load.config('tests/fixtures/test-bisseness.yaml')
   #config <- set.model(config, 'bisseness')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'bisseness')
   #expected <- tree.bisseness(c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01, 0.05, 0.1,
  #                              0.05, 0.1), max.taxa=20)
   #checkEquals(expected, result)

   ## Test Birth-Death Model
   #config <- load.config('tests/fixtures/test-bd.yaml')
   #config <- set.model(config, 'bd')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'bd')
   #expected <- tree.bd(c(0.1, 0.003), max.taxa=20)
   #checkEquals(expected, result)

   ## Test ClaSSE Model
   #config <- load.config('tests/fixtures/test-classe.yaml')
   #config <- set.model(config, 'classe')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'classe')
   #expected <- tree.classe(c(2.5, 0.5, 0, 0, 0, 5, 2.41, 5.24, 0.5, 0),
  #                           max.taxa=20)
   #checkEquals(expected, result)

   ## Test GeoSSE Model
   #config <- load.config('tests/fixtures/test-geosse.yaml')
   #config <- set.model(config, 'geosse')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'geosse')
   #expected <- tree.geosse(c(1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), max.taxa=20)
   #checkEquals(expected, result)

   ## Test MuSSE Model
   #config <- load.config('tests/fixtures/test-musse.yaml')
   #config <- set.model(config, 'musse')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'musse')
   #expected <- tree.musse(c(0.1, 0.07, 0.2, 0.03, 0.01, 0.06, 0.01, 0.03, 0.04,
  #                           0.08, 0.1, 0.01), max.taxa=20)
   #checkEquals(expected, result)

   #TODO: Complete YAML file for QuaSSE (figure out prior for drift/diffusion)
   ## Test QuaSSE Model
   #config <- load.config('tests/fixtures/test-quasse.yaml')
   #config <- set.model(config, 'quasse')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'quasse')
   #lambda <- function(x) constant.x(x, 0.1)
   #mu <- function(x) constant.x(x, 0.03)
   #char <- make.brownian.with.drift(0, 0.025)
   #expected <- tree.quasse(c(lambda, mu, char), max.taxa=20)
   #checkEquals(expected, result)

   ## Test Yule Models
   #config <- load.config('tests/fixtures/test-yule.yaml')
   #config <- set.model(config, 'yule')
   #theta <- sample.priors(config)
   #result <- speciation.model(theta, nsim=1, tips=20, 'yule')
   #expected <- tree.yule(c(0.1), max.taxa=20)
   #checkEquals(expected, result)

   # Test Multiple Simulations
   config <- load.config('tests/fixtures/test-bisse.yaml')
   config <- set.model(config, 'bisse')
   theta <- sample.priors(config)
   trees <- speciation.model(theta, nsim=50, tips=20, 'bisse')
   result <- length(trees) # 50 trees
   expected <- 50
   checkEquals(expected, result)

   ## Test varying tip numbers
   trees <- speciation.model(theta, nsim=1, tips=5, 'bisse')
   result <- length(trees$tip.label) # 5 tips
   expected <- 5
   checkEquals(expected, result)

   trees <- speciation.model(theta, nsim=1, tips=100, 'bisse')
   result <- length(trees$tip.label) # 100 tips
   expected <- 100
   checkEquals(expected, result)
}
