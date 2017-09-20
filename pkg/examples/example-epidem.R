require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-epidem.yaml')
config <- set.model(config, 'epidemic')

# simulate target tree
theta <- c(t.end=0.2, N=10000, beta=0.0135, gamma=5, phi=5)
set.seed(50)
obs.tree <- epidem.model(theta, nsim=1, tips=100, model='epidemic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', nthreads=1, model="epidemic", seed=NA, verbose=TRUE)   
