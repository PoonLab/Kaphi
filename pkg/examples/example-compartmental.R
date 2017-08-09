require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

# simulate target tree
theta <- c(t.end=200, N=10000, beta=0.1, gamma=0.002, mu=0.0001, alpha=5)
set.seed(40)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic', fgyResolution=1000)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', model="sir.dynamic", seed=NA, verbose=TRUE)   


# let's examine the contents of the trace file
trace <- read.table('pkg/examples/example-compartmental.tsv', header=T, sep='\t')

# trajectory of mean estimate of beta
par(mar=c(5,5,2,2))
plot(
  sapply(split(trace$beta*trace$weight, trace$n), sum), 
  type='o',
  xlab='Iteration', 
  ylab='Mean beta',
  cex.lab=1,
  main='Trajectory of Mean Beta (SIRD Model, 500 particles)'
  #,ylim=c(0.05,0.105)
)
abline(h=0.1, lty=2)


# use kernel densities to visualize posterior approximations
pal <- rainbow(n=8, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
plot(density
     (trace$beta[trace$n==1], 
     weights=trace$weight[trace$n==1]), 
     col=pal[1], 
     lwd=2, 
     main='SIRD (gamma: shape=1, rate=5)', 
     xlab='SIRD rate parameter (beta)', 
     cex.lab=1.2
     )

for (i in 1:7) {
  temp <- trace[trace$n==i*10,]
  lines(density(temp$beta, weights=temp$weight), col=pal[i+1], lwd=1.5)
}
lines(density
      (trace$beta[trace$n==max(trace$n)], 
      weights=trace$weight[trace$n==max(trace$n)]), 
      col='black', 
      lwd=2
      )  # final estimates
abline(v=0.1, lty=3, col='red')


# show the prior distribution
x <- seq(0, 2, 0.01)
y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[["beta"]]))}
lines(x, y(x), lty=5)

# show posterior distribution (work in progress)
node.heights <- rev(branching.times(obs.tree))

# make a legend
legend(x=1, y=10, legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=60', 'n=70','n=71(final)', 'true beta(0.1)'), lty=c(5,rep(1,9),3), col=c('black', pal, 'black', 'red'), lwd=c(1,2,rep(1.5,7),2,0.75), seg.len=2)



# trajectory of mean estimate of gamma
par(mar=c(5,5,2,2))
plot(
  sapply(split(trace$gamma*trace$weight, trace$n), sum), 
  type='o',
  xlab='Iteration', 
  ylab='Mean gamma',
  cex.lab=1,
  main='Trajectory of Mean Gamma (SIRD Model, 20 particles)'
)
abline(h=1/520, lty=2)

# trajectory of mean estimate of mu
par(mar=c(5,5,2,2))
plot(
  sapply(split(trace$mu*trace$weight, trace$n), sum), 
  type='o',
  xlab='Iteration', 
  ylab='Mean mu',
  cex.lab=1,
  main='Trajectory of Mean Mu (SIRD Model, 20 particles)'
)
abline(h=1/3640, lty=2)



