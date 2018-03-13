require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
#config <- load.config('pkg/examples/example-bisse.yaml')
config <- load.config('tests/fixtures/test-bisse.yaml')
config <- set.model(config, 'bisse')

# simulate target tree
theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=500, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# Visualize parameter identifiability

# calculate kernel distances for varying q01 and q10
png(filename = "~/Documents/BiSSE/corr.kernel500.1.png",width=900,height=900,res=120)
x <- runif(1000, 0, 1)
y <- runif(1000, 0, 1)
res <- mapply(function(val1,val2) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=val1, q10=val2)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val1,val2, "\n")
  mean(dists)
},x,y)
cex.val <- 10*res
# generate a plot
par(mar=c(5,5,2,2))
plot(x, y, type='p', xlab='q01', ylab='q10', cex = cex.val, 
     main='Correlation of q01 and q10')

dev.off()