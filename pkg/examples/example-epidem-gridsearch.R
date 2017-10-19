#------------------------------------------------------------------------------------------------------------
# Grid search for varying parameters in parameter space of epidem.model$theta
require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-epidem.yaml')
config <- set.model(config, 'epidemic')

# simulate target tree
theta <- c(t.end=0.2, N=10000, beta=1, gamma=5, phi=5)
#set.seed(50)
obs.tree <- epidem.model(theta, nsim=1, tips=100, model='epidemic', tsample=0.1, seed=50)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

t.end <- seq(0.05, 0.4, 0.05)
N <- seq(7000, 14000, 1000)
beta <- seq(0.4, 1.8, 0.2)
gamma <- seq(4.4, 5.8, 0.2)
phi <- seq(4.4, 5.8, 0.2)

grid.params <- list('t.end'=t.end, 'N'=N, 'beta'=beta, 'gamma'=gamma, 'phi'=phi)

# create 5-D array and label row names and col names
grid <- array(dim=c(length(t.end), length(N), length(beta), length(gamma), length(phi)), dimnames=list(t.end,N,beta,gamma,phi))
all.combns <- expand.grid(grid.params)

for (param1 in names(grid.params)) {
  for (param2 in names(grid.params)) {
    if (param1 != param2) {
      for (i in grid.params[[param1]]) {
      
        res <- sapply(grid.params[[param2]], function(j) {
          new.theta <- .modify.theta(theta, param1, i, param2, j)
          sim.trees <- epidem.model(new.theta, nsim=5, tips=100, model='epidemic')
          dists <- sapply(sim.trees, function(st) {
            pt <- .preprocess.tree(st, config)
            distance(obs.tree, pt, config)
            
          if (param1 == 't.end' && param2 == 'N') {
            grid[i, j, ]
          }
          
          })
          
          
        })
      }
    }
  }
}


.modify.theta <- function(theta, param1, i, param2, j) {
  new.theta <-sapply(names(theta), function(x) {
    if (param1 == x) {
      theta[[x]] <- i
    } else if (param2 == x) {
      theta[[x]] <- j
    } else {
      theta[[x]] <- theta[[x]]
    }
  })
  return(new.theta)
}