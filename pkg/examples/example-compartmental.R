require(Kaphi)

setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.nondynamic')

config$nsample <- 5

# simulate target tree
#theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
theta <- c(t.end=1000, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
set.seed(7)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.nondynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# calculate kernel distances for varying t.end
x <- seq(100, 1300, 500)    # (from, to, step)
res <- sapply(x, function(value) {
  theta <- c(t.end=value, N=1000, beta=0.01, gamma=1/520, mu=1/3640, alpha=0)
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
plot(x, res, type='b', xlab='t.end', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,1))



## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)

# this takes about....idk how long to run
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', model="sir.nondynamic")    #require a tsv file here ... find a dataset for this epidemiological model

# examine the contents of the trace file
trace <- read.table('pkg/examples/example-compartmental.tsv', header=T, sep='\t')

# trajectory of mean estimate of ...?