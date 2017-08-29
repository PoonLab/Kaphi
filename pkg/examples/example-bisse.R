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

## now let's estimate that posterior distribution!
# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-bisse2.tsv', model='bisse', verbose=TRUE)

# let's examine the contents of the trace file
trace <- read.table('pkg/examples/example-bisse2.tsv', header=T, sep='\t')

#------------------------------------------------------------------------------
# Plot trajectory of mean estimate of lambda and mu

pal <- rainbow(n=6, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
#png('bd-mean01.png')

plot(
  sapply(split(trace$lambda0*trace$weight, trace$n), sum), 
  ylim=c(0, 1), 
  type='l',
  xlab='Iteration', 
  ylab='Mean Parameter Value',
  cex.lab=1,
  main='Trajectory of Mean Lambda and Mu (BiSSE Model, 100 particles)',
  col=pal[1]
)
lines(
  sapply(split(trace$mu0*trace$weight, trace$n), sum),
  type='l',
  col=pal[2]
)
lines(
  sapply(split(trace$lambda1*trace$weight, trace$n), sum),
  type='l',
  col=pal[3]
)
lines(
  sapply(split(trace$mu1*trace$weight, trace$n), sum),
  type='l',
  col=pal[4]
)
lines(
  sapply(split(trace$q01*trace$weight, trace$n), sum),
  type='l',
  col=pal[5]
)
lines(
  sapply(split(trace$q10*trace$weight, trace$n), sum),
  type='l',
  col=pal[6]
)
abline(h=0.1, lty=2, col=pal[1])
abline(h=0.003, lty=2, col=pal[2])
abline(h=0.1, lty=2, col=pal[3])
abline(h=0.003, lty=2, col=pal[4])
abline(h=0.01, lty=2, col=pal[5])
abline(h=0.01, lty=2, col=pal[6])
#dev.off()


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
