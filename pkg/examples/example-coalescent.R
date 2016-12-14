require(Kaphi)

# load configuration file
config <- load.config('pkg/examples/example-coalescent.yaml')
config <- set.model(config, const.coalescent)

# simulate target tree
theta <- c(Ne.tau=1000)
set.seed(100)
obs.tree <- const.coalescent(theta, nsim=1, n.tips=20)[[1]]
