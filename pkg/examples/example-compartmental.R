require(Kaphi)

setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-compartmental.yaml')
config <- set.model(config, 'sir.dynamic')

config$nsample <- 5

#simulate target tree
theta <- c(t.end=30.*52, N=1000, beta=0.01, gamma=1/520, mu=1/3640, epsilon=0)
set.seed(7)
obs.tree <- compartmental.model(theta, nsim=1, tips=100, model='sir.dynamic')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)
