require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

# simulate target tree
theta <- c(t.end=200, N=10000, beta=0.1, gamma=0.002, mu=0.0001, alpha=5)
set.seed(50)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic', fgyResolution=1000)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)
result <- run.smc(ws, trace.file='pkg/examples/example-compartmental.tsv', nthreads=1, model="sir.dynamic", seed=NA, verbose=TRUE)   

########################################################################################################################################################################
