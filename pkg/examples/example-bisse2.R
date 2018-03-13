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
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

# Visualize parameter identifiability

# calculate kernel distances for varying lambda0
png(filename = "~/Documents/BiSSE/kernel50_1_Lambda0.png",width=900,height=900,res=120)
x <- runif(1000, 0, 1)
res <- sapply(x, function(val) {
  theta <- c(lambda0=val, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='Lambda0', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of Lambda0 (BiSSE Model)')
abline(v=0.2, lty=2)

dev.off()

# calculate kernel distances for varying lambda1
png(filename = "~/Documents/BiSSE/kernel50_1_Lambda1.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.5)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=val, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='Lambda1', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of Lambda1 (BiSSE Model)')
abline(v=0.1, lty=2)

dev.off()

# calculate kernel distances for varying mu0
png(filename = "~/Documents/BiSSE/kernel50_1_mu0.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=val, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='mu0', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of mu0 (BiSSE Model)')
abline(v=0.003, lty=2)

dev.off()

# calculate kernel distances for varying mu1
png(filename = "~/Documents/BiSSE/kernel50_1_mu1.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=val, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='mu1', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of mu1 (BiSSE Model)')
abline(v=0.003, lty=2)

dev.off()

# calculate kernel distances for varying q01
png(filename = "~/Documents/BiSSE/kernel50_1_q01.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.1)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=val, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='q01', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of q01 (BiSSE Model)')
abline(v=0.03, lty=2)

dev.off()

# calculate kernel distances for varying q10
png(filename = "~/Documents/BiSSE/kernel50_1_q10.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.05)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=val)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
  dists <- sapply(sim.trees, function(st) {
    pt <- .preprocess.tree(st, config)
    distance(obs.tree, pt, config)
  })
  cat(val, "\n")
  mean(dists)
})
o <- order(x)
# generate a plot
par(mar=c(5,5,2,2))
plot(x[o], res[o], type='p', xlab='q10', ylab='Mean distance', cex.lab=1.2, ylim=c(0,0.3), 
     main='Identifiability of q10 (BiSSE Model)')
abline(v=0.01, lty=2)

dev.off()