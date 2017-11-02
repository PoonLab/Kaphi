require(Kaphi)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-epidem.yaml')
config <- set.model(config, 'epidemic')

# simulate target tree
theta <- c(t.end=20, N=1000, beta=0.001, gamma=0.3, phi=0.15)
#set.seed(50)
obs.tree <- epidem.model(theta, nsim=1, tips=100, model='epidemic', seed=50)[[1]]
obs.tree <- parse.input.tree(obs.tree, config)


## estimate posterior distribution

# initialize workspace
ws <- init.workspace(obs.tree, config)
result <- run.smc(ws, trace.file='pkg/examples/example-epidem.tsv', nthreads=1, model="epidemic", seed=NA, verbose=TRUE)   



trace <- read.table('pkg/examples/example-epidem.tsv', header=T, sep='\t')

for (param in names(theta)) {
  par(mar=c(5,5,2,2))
  plot(
    sapply(split(trace[[param]]*trace$weight, trace$n), sum), 
    type='o',
    xlab='Iteration', 
    ylab=paste0('Mean ', param),
    cex.lab=1,
    main=paste0('Trajectory of Mean ', param, ' (MASTER, ', config$nparticle, ' particles)')
    #,ylim=c(0.05,0.105)
  )
  # true param value
  abline(h=theta[[param]], lty=2)
  
  
  # use kernel densities to visualize posterior approximations
  pal <- rainbow(n=8, start=0, end=0.5, v=1, s=1)
  par(mar=c(5,5,2,2))
  plot(density
       (trace[[param]][trace$n==1], 
         weights=trace$weight[trace$n==1]), 
       col=pal[1], 
       lwd=2, 
       main=paste0('MASTER ', config$priors[[param]]), 
       xlab=paste0('MASTER parameter (', param, ')'), 
       cex.lab=1.2
  )
  
  for (i in 1: ( length(unique(trace$n)) %/% 10 ) ) {
    temp <- trace[trace$n==i*10,]
    lines(density(temp[[param]], weights=temp$weight), col=pal[i+1], lwd=1.5)
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
  
  # show posterior distribution (work in progress)
  node.heights <- rev(branching.times(obs.tree))
  
  # make a legend
  legend(x=1, y=10, legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=60', 'n=70','n=71(final)', paste0('true ', param, '(', theta[[param]], ')')), lty=c(5,rep(1,9),3), col=c('black', pal, 'black', 'red'), lwd=c(1,2,rep(1.5,7),2,0.75), seg.len=2)
}

