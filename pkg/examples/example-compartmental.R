require(Kaphi)

setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')
#config$nsample <- 1
# simulate target tree
theta <- c(t.end=200, N=5000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
set.seed(40)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic', fgyResolution=1000)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)

result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', model="sir.dynamic", seed=NA, verbose=TRUE)   

# let's examine the contents of the trace file
trace <- read.table('pkg/examples/example-compartmental.tsv', header=T, sep='\t')

# trajectory of mean estimate of lambda
par(mar=c(5,5,2,2))
#png('compartmental-2np.png')
plot(
  sapply(split(trace$beta*trace$weight, trace$n), sum), 
  type='o',
  xlab='Iteration', 
  ylab='Mean beta',
  cex.lab=1,
  main='Trajectory of Mean Beta (SIRD Model, 2 particles)'
)
abline(h=0.1, lty=2)
#abline(h=0.09, lty=2)
#abline(h=0.11, lty=2)
#dev.off()

# use kernel densities to visualize posterior approximations
pal <- rainbow(n=8, start=0, end=0.5, v=1, s=1)
par(mar=c(5,5,2,2))
png('yule-1000-dens_3.png')
plot(density(trace$lambda[trace$n==1], weights=trace$weight[trace$n==1]), xlim=c(0, 2), col=pal[1], lwd=2, main='Yule (gamma: shape=2, rate=1)', xlab='Yule rate parameter (lambda)', cex.lab=1.2, ylim=c(0, 15))

for (i in 1:7) {
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
legend(x=1, y=10, legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=60', 'n=70','n=71(final)', 'true lambda(0.1)'), lty=c(5,rep(1,9),3), col=c('black', pal, 'black', 'red'), lwd=c(1,2,rep(1.5,7),2,0.75), seg.len=2)
dev.off()



#########################################################################
# calculate kernel distances for varying N
x <- seq(1000,10000, 1000)    # (from, to, step)
resx <- sapply(x, function(value) {
  theta <- c(t.end=50, N=value, beta=0.1, gamma=1/520., mu=1/3640., alpha=0)
  sim.trees <- compartmental.model(theta, nsim=100, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, resx, type='b', xlab='N', ylab='Mean kernel distance', cex.lab=1.2)


# calculate kernel distances for varying beta
y <- seq(0.01,0.235, 0.025)    # (from, to, step)
resy <- sapply(y, function(value) {
  theta <- c(t.end=50, N=5000, beta=value, gamma=1/520, mu=1/3640, alpha=0)
  sim.trees <- compartmental.model(theta, nsim=100, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(y, resy, type='b', xlab='beta', ylab='Mean kernel distance', cex.lab=1.2)


#calculate kernel distances for varying gamma
x <- seq(0.0005,0.03, 0.0015)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=50, N=5000, beta=0.1, gamma=value, mu=1/3640, alpha=0)
  sim.trees <- compartmental.model(theta, nsim=100, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='b', xlab='gamma', ylab='Mean kernel distance', cex.lab=1.2)


#calculate kernel distances for varying mu
x <- seq(0.00005,0.0006, 0.00005)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=value, alpha=0)
  sim.trees <- compartmental.model(theta, nsim=100, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='b', xlab='mu', ylab='Mean kernel distance', cex.lab=1.2)


#calculate kernel distances for varying alpha
#change test to seir (uses alpha parameter) and resimulate target tree
config <- set.model(config, 'seir')
theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=5)
set.seed(50)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='seir')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

x <- seq(0.0001,0.01, 0.0005)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=value)
  sim.trees <- compartmental.model(theta, nsim=10, tips=100, model='seir')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})
# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='b', xlab='alpha', ylab='Mean kernel distance', cex.lab=1.2)



#--------------------------------------------------------------------
# Grid search for all pairwise combinations of values {beta} x {N}

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=0)
set.seed(25)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

x <- seq(1000,10000, 1000)   # N
y <- seq(0.01,0.235, 0.025)  # beta

m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x,y))
ind <- 1

# fill columns
for (i in y) {
  cat('beta: ', i, '\n')
  res <- sapply(x, function(val) {
    theta <- c(t.end=50, N=val, beta=i, gamma=1/520, mu=1/3640, alpha=0)
    sim.trees <- compartmental.model(theta, nsim=50, tips=50, model='sir.dynamic')
    dists <- sapply(sim.trees, function(st) {
      pt <- .preprocess.tree(st, config)
      distance(obs.tree, pt, config)
    })
    cat(' N: ', val, "\n")
    mean(dists)
  })
  cat('writing values to col. ', ind, '\n')
  m[,ind] <- res
  ind <- ind + 1
}

# plot heatmap using gplots
require(grDevices)
require(gplots)
pal <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
hm1 <- heatmap.2(m, 
                Rowv=NA, Colv=NA,
                scale="none", na.rm=TRUE,
                col=pal,  
                margins=c(5,5),
                trace='none',
                density.info='none',
                ylab='N', xlab='beta', 
                main='Distance from obs.tree')

#--------------------------------------------------------------------
# Grid search for all pairwise combinations of values {gamma} x {mu}

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=0)
set.seed(25)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

x <- seq(0.0005,0.0075, 0.00075)   # gamma
y <- seq(0.00005,0.00055, 0.00005)  # mu

m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x,y))
ind <- 1

# fill columns
for (i in y) {
  cat('mu: ', i, '\n')
  res <- sapply(x, function(val) {
    theta <- c(t.end=50, N=5000, beta=0.1, gamma=val, mu=i, alpha=0)
    sim.trees <- compartmental.model(theta, nsim=50, tips=50, model='sir.dynamic')
    dists <- sapply(sim.trees, function(st) {
      pt <- .preprocess.tree(st, config)
      distance(obs.tree, pt, config)
    })
    cat(' gamma: ', val, "\n")
    mean(dists)
  })
  cat('writing values to col. ', ind, '\n')
  m[,ind] <- res
  ind <- ind + 1
}

# plot heatmap using gplots
require(grDevices)
require(gplots)
pal <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
hm1 <- heatmap.2(m, 
                 Rowv=NA, Colv=NA,
                 scale="none", na.rm=TRUE,
                 col=pal,  
                 margins=c(5,5),
                 trace='none',
                 density.info='none',
                 ylab='gamma', xlab='mu', 
                 main='Distance from obs.tree')