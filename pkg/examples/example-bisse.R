require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
#config <- load.config('pkg/examples/example-bisse.yaml')
config <- load.config('tests/fixtures/test-bisse.yaml')
config <- set.model(config, 'bisse')

# simulate target tree
theta <- c(lambda0=0.1, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.01, q10=0.01)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# calculate kernel distances for varying lambda0
x <- seq(0.01, 0.3, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda0=val, lambda1=0.1, mu0=0.003, mu1=0.03, q01=0.01, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=50, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='o', xlab='Lambda0', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of Lambda0 (BiSSE Model)')
abline(v=0.1, lty=2)
