## Measure the kernel distance function for a Yule tree
require(Kaphi)
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
config <- load.config('pkg/examples/example-yule.yaml')
config <- set.model(config, 'yule')

# simulate target tree
theta <- c(lambda=0.1)  # this is the true value
set.seed(51)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='yule')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!

# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-yule2.tsv', nthreads = 10, model='yule', verbose=TRUE)

# let's examine the contents of the trace file

trace <- read.table('pkg/examples/example-yule2.tsv', header=T, sep='\t')

pdf(file='~/Documents/exyule.kernel.1.pdf')

# trajectory of mean estimate of lambda
for (param in names(theta)) {
  par(mar=c(5,5,2,2))
  plot(
    sapply(split(trace[[param]]*trace$weight, trace$n), sum), 
    ylim=c(0, 2), 
    type='o',
    xlab='Iteration', 
    ylab=paste0('Mean ', param),
    cex.lab=1,
    main=paste0('Trajectory of Mean ', param, ' (Yule Model, ', config$nparticle, ' particles)')
  )
  #true param value
  abline(h=theta[[param]], lty=2)
  
  # use kernel densities to visualize posterior approximations
  pal <- rainbow(n=(length(unique(trace$n)) %/% 10)+1, start=0, end=0.5, v=1, s=1)
  par(mar=c(5,5,2,2))
  plot(density
       (trace[[param]][trace$n==1], 
         weights=trace$weight[trace$n==1]), 
       col=pal[1], 
       lwd=2, 
       main=paste0('Yule ', config$priors[[param]]), 
       xlab=paste0('Yule rate parameter (', param, ')'), 
       cex.lab=1.2
  )
  
  for (i in ( length(unique(trace$n)) %/% 10 ) ) {
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
  legend(x=1, y=10, legend=c('prior', 'n=1', 'n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'n=53(final)', 'true lambda(0.1)'), lty=c(5,rep(1,7),3), col=c('black', pal, 'black', 'red'), lwd=c(1,2,rep(1.5,4),2,0.75), seg.len=2)
}

dev.off()
