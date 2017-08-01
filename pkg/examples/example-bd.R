## Example of Kaphi with Birth-Death speciation model
require(Kaphi)

# if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-bd.yaml')
config <- set.model(config, 'bd')

# simulate target tree
theta <- c(lambda=0.1, mu=0.003)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bd')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!
# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-bd2.tsv', model='bd', verbose=TRUE)

# let's examine the contents of the trace file
trace <- read.table('pkg/examples/example-bd2.tsv', header=T, sep='\t')


#------------------------------------------------------------------------------
# Plot trajectory of mean estimate of lambda and mu

par(mar=c(5,5,2,2))
#png('bd-mean01.png')
plot(
  sapply(split(trace$lambda*trace$weight, trace$n), sum), 
  ylim=c(0, 1), 
  type='o',
  xlab='Iteration', 
  ylab='Mean Parameter Value',
  cex.lab=1,
  main='Trajectory of Mean Lambda and Mu (Birth-Death Model, 10 particles)'
)
lines(
  sapply(split(trace$mu*trace$weight, trace$n), sum),
  type='o',
  col='red'
)
abline(h=0.1, lty=2, col='black')
abline(h=0.003, lty=2, col='red')
#dev.off()


#------------------------------------------------------------------------------
# use kernel densities to visualize posterior approximations
pal <- rainbow(n=9, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
#png('bd-dist01.png')
plot(
  density(trace$lambda[trace$n==1], weights=trace$weight[trace$n==1]), 
  xlim=c(0, 2), 
  col=pal[1], 
  lwd=2, 
  main='Birth-Death', 
  xlab='Lambda', 
  cex.lab=1.2, 
  ylim=c(0, 15)
)

for (i in 1:8) {
  temp <- trace[trace$n==i*10,]
  lines(density(temp$lambda, weights=temp$weight), 
        col=pal[i+1], lwd=1.5)
}
# final estimates
lines(density(trace$lambda[trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), 
      col='black', lwd=2)  
abline(v=0.1, lty=3, col='red')

# show the prior distribution
x <- seq(0, 2, 0.01)
y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[["lambda"]]))}
lines(x, y(x), lty=5)

# show posterior distribution (work in progress)
node.heights <- rev(branching.times(obs.tree))

# make a legend
legend(
  x=1, y=15, 
  legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=60', 'n=70', 'n=80', 'n=89(final)', 'true lambda(0.1)'), 
  lty=c(5,rep(1,10),3), 
  col=c('black', pal, 'black', 'red'), 
  lwd=c(1,2,rep(1.5,6),2,0.75), 
  seg.len=2
)
#dev.off()


#------------------------------------------------------------------------------
# use kernel densities to visualize posterior approximations
pal <- rainbow(n=9, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
#png('bd-dist01.png')
plot(
  density(trace$mu[trace$n==1], weights=trace$weight[trace$n==1]), 
  xlim=c(0, 2), 
  col=pal[1], 
  lwd=2, 
  main='Birth-Death', 
  xlab='Mu', 
  cex.lab=1.2, 
  ylim=c(0, 15)
)

for (i in 1:8) {
  temp <- trace[trace$n==i*10,]
  lines(density(temp$mu, weights=temp$weight), 
        col=pal[i+1], lwd=1.5)
}
# final estimates
lines(density(trace$mu[trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), 
      col='black', lwd=2)  
abline(v=0.003, lty=3, col='red')

# show the prior distribution
x <- seq(0, 2, 0.01)
y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[["mu"]]))}
lines(x, y(x), lty=5)

# show posterior distribution (work in progress)
node.heights <- rev(branching.times(obs.tree))

# make a legend
legend(
  x=1, y=15, 
  legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=60', 'n=70', 'n=80', 'n=89(final)', 'true mu(0.003)'), 
  lty=c(5,rep(1,10),3), 
  col=c('black', pal, 'black', 'red'), 
  lwd=c(1,2,rep(1.5,6),2,0.75), 
  seg.len=2
)
#dev.off()


#------------------------------------------------------------------------------
# Visualize parameter identifiability

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
y <- seq(0, 0.05, 0.001)
res2 <- sapply(y, function(val) {
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
plot(y, res2, type='o', xlab='Mu', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0.05,0.09),
     main='Identifiability of Mu (Birth-Death Model)')
abline(v=0.003, lty=2)
# log transformed plot
plot(log(y), res2, type='o', xlab='Mu', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0.05,0.12),
     main='Identifiability of Mu (Birth-Death Model)')
abline(v=log(0.003), lty=2)


#------------------------------------------------------------------------------
# Grid search + heatmap for all pairwise combinations of values {lambda} x {mu}

# set up matrix
x <- seq(0.05, 0.3, 0.005)
y <- seq(0, 0.05, 0.001)
m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x,y))
ind <- 1

# fill matrix with distances from obs.tree for each pairwise combination 
#   of lambda and mu.
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

# plot heat map using gplots
require(grDevices)
require(gplots)
pal <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
hm1 <- heatmap.2(m, 
                Rowv=NA, Colv=NA,
                scale="none", na.rm=TRUE,
                col=pal,  
                margins=c(5,5),
                density.info='none',
                ylab='Lambda', xlab='Mu', 
                main='Distance from obs.tree',
                trace='none',
                add.expr=abline(v=c(3.75), h=c(41), lty=2)
                )
