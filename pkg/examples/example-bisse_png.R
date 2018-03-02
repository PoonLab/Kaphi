require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
#config <- load.config('pkg/examples/example-bisse.yaml')
config <- load.config('tests/fixtures/test-bisse (lognorm) unlabelled.yaml')
config <- set.model(config, 'bisse')

# simulate target tree
theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=500, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!
# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-bisse2.tsv', model='bisse', nthreads=10, verbose=TRUE)

#------------------------------------------------------------------------------
trace <- read.table('pkg/examples/example-bisse2.tsv', header=T, sep='\t')

for (param in names(theta)) {
  png(filename = paste0("~/Documents/BiSSE/unlabelled_500_1_", param, ".png"),width=900,height=900,res=120)
  # use kernel densities to visualize posterior approximations
  pal <- rainbow(n=(length(unique(trace$n)) %/% 10)+1, start=0, end=0.5, v=1, s=1)
  plot(density
       (trace[[param]][trace$n==1], 
         weights=trace$weight[trace$n==1]), 
       col=pal[1], 
       lwd=2, 
       main=paste0(' (BiSSE ', config$priors[[param]]), 
       xlab=paste0(' (BiSSE rate parameter (', param, ')',
                   '\nMean: ',
                   mean(trace[[param]][trace$n==max(trace$n)]), 
                   '    Median: ', 
                   median(trace[[param]][trace$n==max(trace$n)]),
                   '\n95% CI (',
                   quantile(trace[[param]][trace$n==max(trace$n)], c(0.025, 0.975))[1],
                   ' , ',
                   quantile(trace[[param]][trace$n==max(trace$n)], c(0.025, 0.975))[2],
                   ')'), 
       cex.lab=0.8
  )
  
  if ( length(unique(trace$n)) >= 10) {
    for (i in 1: ( length(unique(trace$n)) %/% 10 ) ) {
      temp <- trace[trace$n==i*10,]
      lines(density(temp[[param]], weights=temp$weight), col=pal[i+1], lwd=1.5)
    }
  }
  lines(density
        (trace[[param]][trace$n==max(trace$n)], 
          weights=trace$weight[trace$n==max(trace$n)]), 
        col='black', 
        lwd=2
  )  # final estimates
  abline(v=theta[[param]], lty=3, col='red')
  
  
  # show the prior distribution
  x <- sort( replicate(1000, eval(parse(text=config$priors[[param]]))) )
  y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[[param]]))}
  lines(x, y(x), lty=5)
  
}

dev.off()


# Visualize parameter identifiability

# calculate kernel distances for varying lambda0
png(filename = "~/Documents/BiSSE/500_1_Lambda0.png",width=900,height=900,res=120)
x <- runif(1000, 0, 1)
res <- sapply(x, function(val) {
  theta <- c(lambda0=val, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='Lambda0', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of Lambda0 (BiSSE Model)')
abline(v=0.2, lty=2)

dev.off()

# calculate kernel distances for varying lambda1
png(filename = "~/Documents/BiSSE/500_1_Lambda1.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.5)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=val, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='Lambda1', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of Lambda1 (BiSSE Model)')
abline(v=0.1, lty=2)

dev.off()

# calculate kernel distances for varying mu0
png(filename = "~/Documents/BiSSE/500_1_mu0.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=val, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='mu0', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of mu0 (BiSSE Model)')
abline(v=0.003, lty=2)

dev.off()

# calculate kernel distances for varying mu1
png(filename = "~/Documents/BiSSE/500_1_mu1.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.01)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=val, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='mu1', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of mu1 (BiSSE Model)')
abline(v=0.003, lty=2)

dev.off()

# calculate kernel distances for varying q01
png(filename = "~/Documents/BiSSE/500_1_q01.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.1)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=val, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='q01', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of q01 (BiSSE Model)')
abline(v=0.03, lty=2)

dev.off()

# calculate kernel distances for varying q10
png(filename = "~/Documents/BiSSE/500_1_q10.png",width=900,height=900,res=120)
x <- runif(1000, 0, 0.05)
res <- sapply(x, function(val) {
  theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=val)
  sim.trees <- speciation.model(theta, nsim=1, tips=500, model='bisse')
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
plot(x[o], res[o], type='o', xlab='q10', ylab='Mean kernel distance', cex.lab=1.2, ylim=c(0,0.5), 
     main='Identifiability of q10 (BiSSE Model)')
abline(v=0.01, lty=2)

dev.off()