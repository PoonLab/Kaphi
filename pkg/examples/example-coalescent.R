## Measure the kernel distance function for an unscaled coalescent tree
require(Kaphi)

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-coalescent.yaml')
config <- set.model(config, const.coalescent)
config$rbf.variance <- 100

# simulate target tree
theta <- c(Ne.tau=1000)
set.seed(50)
obs.tree <- const.coalescent(theta, nsim=1, n.tips=100)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# calculate kernel distances for varying Ne.tau
x <- seq(100, 2000, 100)
res <- sapply(x, function(val) {
    theta <- c(Ne.tau=val)
    sim.trees <- const.coalescent(theta, nsim=50, n.tips=100)
    dists <- sapply(sim.trees, function(st) {
        pt <- preprocess.tree(st, config)
        distance(obs.tree, pt, config)
    })
    cat(val, "\n")
    mean(dists)
})

# generate a plot
par(mar=c(5,5,2,2))
plot(x,res, type='b', xlab='Ne.tau', ylab='Mean kernel distance', cex.lab=1.2)


## now let's estimate that posterior distribution!

# initialize workspace
ws <- init.workspace(obs.tree, config)
result <- run.smc(ws, trace.file='pkg/examples/example-coalescent.tsv')
