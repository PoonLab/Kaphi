#------------------------------------------------------------------------------------------------------------
# Grid search for varying parameters in parameter space of epidem.model$theta
require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

# simulate target tree
theta <- c(t.end=200, N=10000, beta=0.1, gamma=0.2, mu=0.01, alpha=5)
#set.seed(50)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic', fgyResolution=1000, seed=50)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

t.end <- seq(50, 400, 50)
N <- seq(7000, 14000, 1000)
beta <- seq(0.01, 0.36, 0.05)
gamma <- seq(0.01, 0.36, 0.05)
mu <- seq(0.001, 0.071, 0.01)
alpha <- 5

grid.params <- list('t.end'=t.end, 'N'=N, 'beta'=beta, 'gamma'=gamma, 'mu'=mu, 'alpha'=alpha)

# create 5-D array and label row names and col names
grid <- array(dim=c(length(t.end), length(N), length(beta), length(gamma), length(mu), length(alpha)), dimnames=list(t.end,N,beta,gamma,mu,alpha))
all.combns <- expand.grid(grid.params)

for(x in 1:nrow(all.combns)) {
  theta <- all.combns[x,]
  sim.trees <- compartmental.model(theta, nsim=5, tips=100, model='sir.dynamic', fgyResolution=1000)
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  # populate the index in the matrix
  ind1 <- which(t.end == theta$t.end)
  ind2 <- which(N == theta$N)
  ind3 <- which(beta == theta$beta)
  ind4 <- which(gamma == theta$gamma)
  ind5 <- which(mu == theta$mu)
  ind6 <- 1
  cat("Populating index: [", theta$t.end, ',', theta$N, ',', theta$beta, ',', theta$gamma, ',', theta$mu, ',', theta$alpha, '] with', mean(dists), '\n')
  cat("Indices and values: [ind1=", ind1, 'ind2=', ind2, 'ind3=', ind3, 'ind4=', ind4, 'ind5=', ind5, 'ind6=', ind6, '\n\n')
  grid[ind1, ind2, ind3, ind4, ind5, ind6] <- mean(dists)
  grid
}

require(reshape2)
run.copy <- grid
melted <- melt(run.copy)
write.csv(melted, file='~/Documents/gridsearch.run1.compart.csv')