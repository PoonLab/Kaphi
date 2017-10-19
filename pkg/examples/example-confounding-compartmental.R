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