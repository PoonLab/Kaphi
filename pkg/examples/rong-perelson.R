require(Kaphi)


# Model parameters
# default values drawn from Rong and Perelson (2009) PLCB e1000533
parms <- list(
	lambda=1e4,  # growth rate of uninfected cells (per mL per day)
	d.T=0.01,    # death rate of uninfected cells (per day)
	k=2.4e-8,    # rate of infection (mL/day)
	eta=0.01,    # probability of entering latent state
	d.0=0.001,   # death rate of latently-infected cells
	a.L=0.1,     # rate of transition from latently to productively infected cells
	delta=1.0,   # death rate of productively infected cells (per day)
	N=2000,      # number of virions produced by cell death
	c=23.        # clearance rate of free virus (per day)
)


demes <- c(
    "V",  # free virus
    "L",  # latent infected cells
    "Ts"  # active infected cells
)

births <- rbind(
	c('0', '0', '0'),
	c('0', '0', '0'),
	c('parms$N*parms$delta*Ts', '0', '0')  # Ts->V
)
rownames(births) <- colnames(births) <- demes

migrations <- rbind(
	c('0', 'parms$eta*parms$k*T*V', '(1-parms$eta)*parms$k*T*V'),
	c('0', '0', 'parms$a.L * L'),
	c('0', '0', '0')
)
rownames(migrations) <- colnames(migrations) <- demes

deaths <- c(
	'parms$c * V',      # free viruses removed at rate (c)
	'parms$d.0 * L',    # latently-infected cells removed at rate (d.0)
	'parms$delta * Ts'  # active infected cells removed at rate by bursting (delta)
)
names(deaths) <- demes

ndd <- c('parms$lambda - parms$d.T * T - parms$k * V * T')
names(ndd) <- c('T')
expr <- parse.ode(births, deaths, ndd, migrations)

# ----------------------------------------------------

simulate.RP <- function(expr, demes, parms, start.time, end.time, integ.method, fgy.resol, x0, nsamples.V, nsamples.T, ntimes) {	
	sol <- solve.ode(expr, t0=start.time, t1=end.time, x0=x0, parms=parms, time.pts=fgy.resol, integrationMethod=integ.method)
	
	# initialize sample state matrix
	sample.states <- matrix(0, ncol=3, nrow=ntimes * nsamples)
	colnames(sample.states) <- demes
	
	# assume uniform sampling over times
	time.points <- seq(fgy.resol, 0, -ceiling(fgy.resol/ntimes))[1:ntimes]
	nsamples <- nsamples.V + nsamples.T
	sample.times <- sol$sol[rep(time.points, each=nsamples), 1]
	
	for (i in 1:length(time.points)) {
		tp <- time.points[i]
		row <- as.list(sol$sol[tp,])  # extract time slice
		nsamples.L <- rbinom(1, nsamples.T, row$L/(row$L+row$Ts))  # number of latent cells in sample as binomial outcome
		nsamples.Ts <- nsamples.T - nsamples.L
		to.col <- c(rep(1, nsamples.V), rep(2, nsamples.L), rep(3, nsamples.Ts))
	
		# use (to.col) to assign 1's to random permutation of rows within time point block
		permut <- sample(1:nsamples, nsamples)
		for (j in 1:nsamples) {
			sample.states[permut[j]+(i-1)*nsamples, to.col[j]] <- 1
		}
	}
	
	tree <- simulate.ode.tree(sol, sample.times, sample.states, simulate.migrations=TRUE)
	
	# apply sample states to label tips
	tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) demes[which(x == 1)]), sep="."); 
	
	# Use sample paths to adjust branch lengths.
	# Assume free virus and virus in latent cells have near-zero evolution 
	#  relative to viruses in actively infected cells.
	tree.2 <- tree
	tree.2$edge.length <- sapply(tree$samplePath, function(x) {
		delta.t <- diff(x[,1])
		states <- (x[,2][1:nrow(x)-1])
		return(sum(delta.t[states==3]))
	})
	return(tree.2)
}



start.time <- 0
end.time <- 300  # time elapsed in units of days

integ.method <- 'rk4'
fgy.resol <- 1e4
x0 <- c(V=50, T=600, Ts=0.3, L=2)

nsamples.V <- 50
nsamples.T <- 50
ntimes <- 1

# batch processing here

