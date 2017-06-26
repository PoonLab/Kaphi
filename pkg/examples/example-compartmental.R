require(Kaphi)

setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.nondynamic')

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=1300, N=500, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
set.seed(50)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.nondynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)



## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)

# this takes about....idk how long to run
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', model="sir.nondynamic")    #require a tsv file here ... find a dataset for this epidemiological model



#########################################################################
# calculate kernel distances for varying N
x <- seq(400, 600, 10)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=1000, N=value, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
  sim.trees <- compartmental.model(theta, nsim=10, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})

# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='b', xlab='N', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0.55,0.65))


# calculate kernel distances for varying beta
x <- seq(0, 4, 0.1)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=1000, N=500, beta=value, gamma=1/520, mu=1/3640, alpha=0)
  sim.trees <- compartmental.model(theta, nsim=10, tips=100, model='sir.nondynamic')
  distances <- sapply(sim.trees, function(singletree) {
    processtree <- .preprocess.tree(singletree, config)
    distance(obs.tree, processtree, config)
  })
  cat(value, "\n")
  mean(distances)
})

# generate a plot
par(mar=c(5,5,2,2))
plot(x, res, type='b', xlab='t.end', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0.55,0.65))
