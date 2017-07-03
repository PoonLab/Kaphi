## Measure the kernel distance function for a Yule tree
require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-yule.yaml')
config <- set.model(config, 'yule')

# simulate target tree
theta <- c(lambda=0.1)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=100, model='yule')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# calculate kernel distances for varying lambda
x <- seq(0.01, 0.3, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda=val)
  sim.trees <- speciation.model(theta, nsim=50, tips=100, model='yule')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})

# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='o', xlab='Lambda', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.2), 
     main='Identifiability of Lambda (Yule Model) - Normalized')
abline(v=0.1, lty=2)

## now let's estimate that posterior distribution!

# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-yule2.tsv', model='yule', verbose=F)


# let's examine the contents of the trace file
trace <- read.table('pkg/examples/example-yule2.tsv', header=T, sep='\t')

# trajectory of mean estimate of lambda
par(mar=c(5,5,2,2))
plot(
  sapply(split(trace$lambda*trace$weight, trace$n), sum), 
  ylim=c(0, 0.8), 
  type='b',
  xlab='Iteration', 
  ylab='Mean lambda',
  cex.lab=1.2
)
abline(h=0.09, lty=2)
abline(h=0.11, lty=2)

# use kernel densities to visualize posterior approximations
pal <- rainbow(n=6, start=0, end=0.3, v=0.8, s=0.5)
par(mar=c(5,5,2,2))
plot(density(trace$lambda[trace$n==1], weights=trace$weight[trace$n==1]), xlim=c(-0.5, 2), col=pal[1], lwd=2, main=NA, xlab='Yule rate parameter (lambda)', cex.lab=1.2, ylim=c(0, 10))
for (i in 1:8) {
  temp <- trace[trace$n==i*2,]
  lines(density(temp$lambda, weights=temp$weight), col=pal[i+1], lwd=1.5)
}
lines(density(trace$lambda[trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), col=pal[length(pal)], lwd=2)  # final estimates

# show the prior distribution
x <- seq(0, 2, 0.01)
y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[["lambda"]]))}
lines(x, y(x), lty=2)

# show posterior distribution (work in progress)
node.heights <- rev(branching.times(obs.tree))

# make a legend
legend(x=2300, y=0.0025, legend=c('prior', 't=1', 't=10', 't=20', 't=30', 't=40', 't=47'), lty=c(2,rep(1,6)), col=c('black', pal), lwd=c(1,2,rep(1.5,4),2), seg.len=3)

