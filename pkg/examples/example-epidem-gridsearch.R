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

t.end <- seq(0.1, 0.8, 0.1)
N <- seq(7000, 14000, 1000)
beta <- seq(0.4, 1.8, 0.2)
gamma <- seq(4.4, 5.8, 0.2)
phi <- seq(4.4, 5.8, 0.2)

grid.params <- list('t.end'=t.end, 'N'=N, 'beta'=beta, 'gamma'=gamma, 'phi'=phi)

# create 5-D array and label row names and col names
grid <- array(dim=c(length(t.end), length(N), length(beta), length(gamma), length(phi)), dimnames=list(t.end,N,beta,gamma,phi))
all.combns <- expand.grid(grid.params)

for(x in 1:nrow(all.combns)) {
  theta <- all.combns[x,]
  sim.trees <- epidem.model(theta, nsim=5, tips=100, model='epidemic', tsample=0.1)
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  # populate the index in the matrix
  ind1 <- which(t.end == theta$t.end)
  ind2 <- which(N == theta$N)
  ind3 <- which(beta == theta$beta)
  ind4 <- which(gamma == theta$gamma)
  ind5 <- which(phi == theta$phi)
  cat("Populating index: [", theta$t.end, ',', theta$N, ',', theta$beta, ',', theta$gamma, ',', theta$phi, '] with', mean(dists), '\n')
  cat("Indices and values: [ind1=", ind1, 'ind2=', ind2, 'ind3=', ind3, 'ind4=', ind4, 'ind5=', ind5, '\n\n')
  grid[ind1, ind2, ind3, ind4, ind5] <- mean(dists)
  grid
}


