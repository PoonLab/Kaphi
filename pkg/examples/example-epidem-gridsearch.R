#------------------------------------------------------------------------------------------------------------
# Grid search for varying parameters in parameter space of epidem.model$theta
require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-epidem.yaml')
config <- set.model(config, 'epidemic')

# simulate target tree
theta <- c(t.end=30, N=1000, beta=0.001, gamma=0.3, phi=0.15)
#set.seed(50)
obs.tree <- epidem.model(theta, nsim=1, tips=100, model='epidemic', seed=50)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

t.end <- seq(20, 55, 5)
N <- seq(100, 10600, 1500)
beta <- seq(0.0005, 0.004, 0.0005)
gamma <- seq(0.1, 0.8, 0.1)
phi <- seq(0.05, 0.4, 0.05)

grid.params <- list('t.end'=t.end, 'N'=N, 'beta'=beta, 'gamma'=gamma, 'phi'=phi)

# create 5-D array and label row names and col names
grid <- array(dim=c(length(t.end), length(N), length(beta), length(gamma), length(phi)), dimnames=list(t.end,N,beta,gamma,phi))
all.combns <- expand.grid(grid.params)

for(x in 1:nrow(all.combns)) {
  theta <- all.combns[x,]
  sim.trees <- epidem.model(theta, nsim=5, tips=100, model='epidemic')
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
write.csv(melted, file='~/Documents/Grid-search/gridsearch.run6.csv')


## code for looking at grid search results from issue #120
gr <- read.csv('~/Documents/Grid-search/gridsearch.run6.csv')
names(gr) <- c('index', 't.end', 'N', 'beta', 'gamma', 'phi', 'distance')

pairs <- combn(names(grid.params), 2)

pdf(file='~/Documents/Grid-search/gridsearch.run6.pdf')
for (i in 1:(length(pairs)/2)) {
  param1 <- pairs[1,i]
  param2 <- pairs[2,i]
  z <- matrix(0, nrow=8, ncol=8)
  for (i in 1:8) {
    for (j in 1:8) {
      x <- unique(gr[[param1]])[i]
      y <- unique(gr[[param2]])[j]
      foo <- gr$distance[gr[[param1]]==x & gr[[param2]]==y]
      z[i,j] <- mean(sort(foo))
    }
  }

  filled.contour(x=unique(gr[[param1]]), y=unique(gr[[param2]]), z, color.palette=terrain.colors,
                 plot.title= title(main=paste0('Contour Plot Distances of Varying ', param1, ' and ', param2), xlab=paste0(param1), ylab=paste0(param2)),
                 plot.axes= {axis(1); axis(2); abline(v=theta[[param1]], h=theta[[param2]])}
                 )
}
dev.off()

