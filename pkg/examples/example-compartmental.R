require(Kaphi)

setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.nondynamic')

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=0)
set.seed(25)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.nondynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)

# this takes about....idk how long to run
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', model="sir.nondynamic")    #require a tsv file here ... find a dataset for this epidemiological model



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
config <- set.model(config, 'sir.nondynamic')

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=50, N=5000, beta=0.1, gamma=1/520, mu=1/3640, alpha=0)
set.seed(25)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.nondynamic')[[1]]
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
    sim.trees <- compartmental.model(theta, nsim=50, tips=50, model='sir.nondynamic')
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