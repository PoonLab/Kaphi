## Measure the kernel distance function for a Yule tree
require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-yule.yaml')
config <- set.model(config, 'yule')

# simulate target tree
theta <- c(lambda=0.1)  # this is the true value
set.seed(51)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='yule')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!

# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-yule2.tsv', nthreads = 1, model='yule', verbose=TRUE)

# let's examine the contents of the trace file

trace <- read.table('pkg/examples/example-yule2.tsv', header=T, sep='\t')

# trajectory of mean estimate of lambda
par(mar=c(5,5,2,2))
#png('yule01.png')
plot(
  sapply(split(trace$lambda*trace$weight, trace$n), sum), 
  ylim=c(0, 2), 
  type='o',
  xlab='Iteration', 
  ylab='Mean lambda',
  cex.lab=1,
  main='Trajectory of Mean Lambda (Yule Model, 1000 particles)'
)
abline(h=0.05, lty=2)
#abline(h=0.09, lty=2)
#abline(h=0.11, lty=2)
#dev.off()

# use kernel densities to visualize posterior approximations
pal <- rainbow(n=6, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
#png('yule01.png')
plot(density(trace$lambda[trace$n==1], weights=trace$weight[trace$n==1]), xlim=c(0, 2), col=pal[1], lwd=2, main='Yule (gamma: shape=2, rate=1)', xlab='Yule rate parameter (lambda)', cex.lab=1.2, ylim=c(0, 15))

for (i in 1:5) {
  temp <- trace[trace$n==i*10,]
  lines(density(temp$lambda, weights=temp$weight), col=pal[i+1], lwd=1.5)
  #cat('iter:', i*10, '\n')
  #Sys.sleep(2)
}
lines(density(trace$lambda[trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), col='black', lwd=2)  # final estimates
abline(v=0.1, lty=3, col='red')

# show the prior distribution
x <- seq(0, 2, 0.01)
y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[["lambda"]]))}
lines(x, y(x), lty=5)

# show posterior distribution (work in progress)
node.heights <- rev(branching.times(obs.tree))

# make a legend
legend(x=1, y=10, legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=53(final)', 'true lambda(0.1)'), lty=c(5,rep(1,7),3), col=c('black', pal, 'black', 'red'), lwd=c(1,2,rep(1.5,4),2,0.75), seg.len=2)
#dev.off()
