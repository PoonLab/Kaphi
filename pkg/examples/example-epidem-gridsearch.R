#------------------------------------------------------------------------------------------------------------
# Grid search for varying parameters in parameter space of epidem.model$theta
require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-epidem.yaml')
config <- set.model(config, 'epidemic')

# simulate target tree
theta <- c(t.end=0.2, N=10000, beta=1, gamma=5, phi=3)
#set.seed(50)
obs.tree <- epidem.model(theta, nsim=1, tips=100, model='epidemic', tsample=0.1, seed=50)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

t.end <- seq(0.1, 1.5, 0.2)
N <- seq(7000, 14000, 1000)
beta <- seq(0.2, 3.2, 0.4)
gamma <- seq(2.8, 5.6, 0.4)
phi <- seq(2.2, 5.0, 0.4)

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


# if R session aborts when attempting to visualize scatter plot, save the data into a file and re-import
require(reshape2)
run3 <- grid
melted <- melt(run3)
write.csv(file='~/Documents/gridsearch.run3.csv')

data <- read.table('~/Documents/gridsearch.run3.csv', sep=',')

t.end <- seq(0.1, 1.5, 0.2)
N <- seq(7000, 14000, 1000)
beta <- seq(0.2, 3.2, 0.4)
gamma <- seq(2.8, 5.6, 0.4)
phi <- seq(2.2, 5.0, 0.4)
grid.params <- list('t.end'=t.end, 'N'=N, 'beta'=beta, 'gamma'=gamma, 'phi'=phi)
grid <- array(dim=c(length(t.end), length(N), length(beta), length(gamma), length(phi)), dimnames=list(t.end,N,beta,gamma,phi))

for (i in 1:length(rownames(data))) {
  params <- data[i,]
  ind1 <- which(as.character(t.end) == as.character(params$Var1))
  ind2 <- which(as.character(N) == as.character(params$Var2))
  ind3 <- which(as.character(beta) == as.character(params$Var3))
  ind4 <- which(as.character(gamma) == as.character(params$Var4))
  ind5 <- which(as.character(phi) == as.character(params$Var5))
  grid[ind1, ind2, ind3, ind4, ind5] <- params$value
}
