require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
#config <- load.config('pkg/examples/example-bisse.yaml')
config <- load.config('tests/fixtures/test-bisse(gamma).yaml')
config <- set.model(config, 'bisse')

# simulate target tree
theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!
# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-bisse2.tsv', model='bisse', nthreads=10, verbose=TRUE)

#------------------------------------------------------------------------------
trace <- read.table('pkg/examples/example-bisse2.tsv', header=T, sep='\t')

pdf(file='~/Documents/exbisse.kernel.1.pdf')           

for (param in names(theta)) {
  par(mar=c(5,5,2,2))
  plot(
    sapply(split(trace[[param]]*trace$weight, trace$n), sum), 
    type='o',
    xlab='Iteration', 
    ylab=paste0('Mean ', param),
    cex.lab=1,
    main=paste0('Trajectory of Mean ', param, ' (BiSSE, ', config$nparticle, ' particles)')
    #,ylim=c(0.05,0.105)
  )
  # true param value
  abline(h=theta[[param]], lty=2)
  
  
  # use kernel densities to visualize posterior approximations
  pal <- rainbow(n=(length(unique(trace$n)) %/% 10)+1, start=0, end=0.5, v=1, s=1)
  par(mar=c(5,5,2,2))
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
