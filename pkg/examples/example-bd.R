require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-bd.yaml')
config <- set.model(config, 'bd')

# simulate target tree
theta <- c(lambda=0.1, mu=0.003)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bd')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# calculate kernel distances for varying lambda
x <- seq(0.01, 0.3, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda=val, mu=0.003)
  sim.trees <- speciation.model(theta, nsim=50, tips=50, model='bd')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='o', xlab='Lambda', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of Lambda (Birth-Death Model)')
abline(v=0.1, lty=2)

# calculate kernel distances for varying mu
y <- seq(0, 0.1, 0.0034)
res <- sapply(y, function(val) {
  theta <- c(lambda=0.1, mu=val)
  sim.trees <- speciation.model(theta, nsim=100, tips=50, model='bd')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(log(y), res, type='o', xlab='Mu', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0.05,0.12),
     main='Identifiability of Mu (Birth-Death Model)')
abline(v=log(0.003), lty=2)

#--------------------------------------------------------------------
# Grid search for all pairwise combinations of values {lambda} x {mu}

# load configuration file
config <- load.config('pkg/examples/example-bd.yaml')
config <- set.model(config, 'bd')

# simulate target tree
theta <- c(lambda=0.1, mu=0.003)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bd')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# set up distance matrix
x <- seq(0.05, 0.3, 0.025)
y <- seq(0, 0.05, 0.005)
m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x,y))
ind <- 1

# fill columns
for (i in y) {
  cat('mu: ', i, '\n')
  res <- sapply(x, function(val) {
    theta <- c(lambda=val, mu=i)
    sim.trees <- speciation.model(theta, nsim=50, tips=50, model='bd')
    dists <- sapply(sim.trees, function(st) {
      pt <- .preprocess.tree(st, config)
      distance(obs.tree, pt, config)
    })
    cat('  lambda: ', val, "\n")
    mean(dists)
  })
  cat('writing values to col. ', ind, '\n')
  m[,ind] <- res
  ind <- ind + 1
}

# plot -- I'll figure out how to do this correctly tomorrow
plot(x, y, cex=sqrt(m))
