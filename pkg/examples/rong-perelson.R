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

simulate.RP <- function(sample.times, is.rna, expr, parms, x0, start.time, end.time, integ.method='rk4', fgy.resol=1e4, max.tries=3) {
	demes <- expr$demeNames
	
	sol <- solve.ode(expr, t0=start.time, t1=end.time, x0=x0, parms=parms, time.pts=fgy.resol, integrationMethod=integ.method)
	
	if (any(is.nan(sol$sol))) {
		# get 10x resolution for first 10%
		t1 <- start.time + 0.1*(end.time-start.time)
		sol0 <- solve.ode(expr, t0=start.time, t1=t1, x0=x0, parms=parms, time.pts=fgy.resol, integrationMethod=integ.method)
		
		# for downsampling 10x interval
		indices <- seq(fgy.resol, 1, -10)
		
		# copy over last entry of 10x interval as steady state solution
		for (i in 1:(fgy.resol-length(indices))) {  # 1..9000
			# i=1 is the last point in time
			sol$F[[i]] <- sol0$F[[1]]
			sol$G[[i]] <- sol0$G[[1]]
			sol$Y[[i]] <- sol0$Y[[1]]
			sol$sol[fgy.resol-(i-1), ] <- sol0$sol[fgy.resol,]
		}
		
		# copy over high-res values
		for (i in length(indices):1) {  # 9001..10000
			sol$F[[fgy.resol-i+1]] <- sol0$F[[indices[i]]]
			sol$G[[fgy.resol-i+1]] <- sol0$G[[indices[i]]]
			sol$Y[[fgy.resol-i+1]] <- sol0$Y[[indices[i]]]
			sol$sol[i,] <- sol0$sol[rev(indices)[i],]
		}
		sol$sol[,1] <- rev(sol$times)
	}
	if (any(is.nan(sol$sol))) {
		# still not a valid ODE solution
		return(NA)
	}
	
	## generate sample states matrix
	
	nsamples <- length(sample.times)
	sample.states <- matrix(0, ncol=3, nrow=nsamples)
	
	samples <- split(1:length(sample.times), sample.times)
	for (i in 1:length(samples)) {
		samp <- samples[[i]]  # vector of tip indices
		
		# locate time point on ODE solution time axis
		time.pt <- as.integer(names(samples)[i])
		h <- as.integer(cut(time.pt, c(-Inf, sol$times)))

		n <- length(samp)
		n.V <- sum(is.rna[samp])
		n.T <- n - n.V
		
		# use numerical solution of ODE to partition cells into L and Ts states
		row <- as.list(sol$sol[h,])  # extract time slice
		n.L <- rbinom(1, n.T, row$L/(row$L+row$Ts))  # number of latent cells in sample as binomial outcome
		n.Ts <- n.T - n.L
		
		col.idx <- c(rep(1, n.V), rep(2, n.L), rep(3, n.Ts))
		col.idx <- sample(col.idx, n)  # permutation
		for (j in 1:n) {
			sample.states[samp[j], col.idx[j]] <- 1
		}
	}
	colnames(sample.states) <- demes
	rownames(sample.states) <- 1:nsamples
	
	# simulate tree from ODE solution
	tries <- 0
	while(tries < max.tries)
	{
		try(tree <- simulate.ode.tree(sol, sample.times, sample.states, simulate.migrations=TRUE))
		 tries <- tries + 1
		if (exists('tree')) {
			break
		}
	}
	if (tries == max.tries) {
		cat("Failed to simulate tree..\n")
		return(NA)
	}
	
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

################################################################

# tree kernel configuration
config <- list(
	decay.factor=0.2, 
	rbf.variance=2.0, 
	sst.control=1, 
	norm.mode='mean'
)

# read tree to compare
setwd('~/Documents/Seminars/Dynamics/Dynamics2017/jones')
obs.tree <- preprocess.tree(read.tree('data/patient_13889.less.rtt.nwk'), config)
obs.labels <- ifelse(grepl("PLASMA", obs.tree$tip.label), 'V', 'C')
obs.denom <- tree.kernel(obs.tree, obs.tree, lambda=config$decay.factor, sigma=config$rbf.variance, label1=obs.labels, label2=obs.labels)


# shift tip dates so that origin roughly coincides with estimated MRCA from 
#  BEAST strict clock analysis - median root height = 6293 days
sample.times <- as.integer(gsub(".+_([0-9]+)", "\\1", obs.tree$tip.label))

# parse tip labels
is.rna <- grepl("PLASMA", obs.tree$tip.label)



# time elapsed in units of days
start.time <- -500  # some tips sampled at t=0
end.time <- max(sample.times)

# initial conditions
x0 <- c(V=1, T=600, Ts=0, L=0)


# ODE settings
fgy.resol <- 1e4  # needs to be cranked up to capture early dynamics
integ.method <- 'rk4'


get.steady.state <- function(params) {
	# Equation (2) from paper (epsilon = 0)
	V.0 <- with(params, N * lambda / c * (1 - d.0 / (d.0 + a.L) * eta) - d.T / k)
	T.0 <- with(params, lambda / (d.T + k * V.0))
	L.0 <- with(params, eta * k * V.0 * T.0 / (d.0 + a.L))
	Ts.0 <- with(params, c * V.0 / (N * delta))
	
	c(V=V.0, T=T.0, L=L.0, Ts=Ts.0)
}

# prior distributions over model parameters
sample.prior <- function(max.tries=1e3) {
	tries <- 0
	while(tries < max.tries) {
		retry <- FALSE
		proposal <- list(
			lambda = rlnorm(1, meanlog=log(1e4), sdlog=3),
			d.T    = rlnorm(1, log(0.01), 3),
			k      = rlnorm(1, log(2.4e-8), 3),
			eta    = rbeta(1, 0.1, 9.9),
			d.0    = rlnorm(1, log(0.001), 3),
			a.L    = rlnorm(1, log(0.1), 3),
			delta  = rlnorm(1, 0, 3),
			N      = rlnorm(1, log(2000), 3),
			c      = rlnorm(1, log(23), 3)
		)
		tries <- tries + 1
		x.inf <- get.steady.state(proposal)
		if (any(x.inf<= 1)) {
			retry <- TRUE
		}
		if (!retry) {
			break
		}
	}
	return (proposal)
}


# prepare output file
cat(c('rep', 'lambda', 'd.T', 'k', 'eta', 'd.0', 'a.L', 'delta', 'N', 'c', 'kernel',  'kernel.norm', 'newick', '\n'), file='trial1.log', sep='\t')





run1 <- function(rep) {
	# sample parameters from prior distribution
	p0 <- sample.prior()

	# simulate tree
	sim.tree <- simulate.RP(sample.times, is.rna, expr, p0, x0, start.time, end.time, integ.method, fgy.resol)
	if (!exists('sim.tree') | is.na(sim.tree)) { return }
	
	# process the resulting tree 
	sim.tree <- preprocess.tree(sim.tree, config)

	sim.labels <- ifelse(grepl("V", sim.tree$tip.label), 'V', 'C')
	sim.denom <- tree.kernel(sim.tree, sim.tree, lambda=config$decay.factor, sigma=config$rbf.variance, label1=sim.labels, label2=sim.labels)

	# compute kernel score
	res <- tree.kernel(obs.tree, sim.tree, lambda=config$decay, sigma=config$rbf, label1=obs.labels, label2=sim.labels)
	res.norm <- res / sqrt(obs.denom * sim.denom)

	# write output line
	newick <- paste0('"', write.tree(sim.tree), '"')
	write(c(rep, p0$lambda, p0$d.T, p0$k, p0$eta, p0$d.0, p0$a.L, p0$delta, p0$N, p0$c, res, res.norm, newick), append=TRUE, file='trial1.log', sep='\t')
}

#require(parallel)
#mclapply(1:1000, function(i) run1(i), mc.cores=4)

for (i in 1:1000) {
	run1(i)
}

