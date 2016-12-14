## Measure the kernel distance function for an unscaled coalescent tree

require(Kaphi)

# load configuration file
config <- load.config('pkg/examples/example-coalescent.yaml')
config <- set.model(config, const.coalescent)
config$rbf.variance <- 100

# simulate target tree
theta <- c(Ne.tau=1000)
set.seed(33)
obs.tree <- const.coalescent(theta, nsim=1, n.tips=100)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

x <- seq(100, 2000, 100)
res <- sapply(x, function(val) {
    theta <- c(Ne.tau=val)
    sim.trees <- const.coalescent(theta, nsim=50, n.tips=100)
    dists <- sapply(sim.trees, function(st) {
        pt <- preprocess.tree(st, config)
        distance(obs.tree, pt, config)
    })
    mean(dists)
})

plot(x, res)
