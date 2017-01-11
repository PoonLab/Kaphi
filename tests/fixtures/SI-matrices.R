# rcolgem definition for a standard SIS model with recovery
# This model has a relatively simple analytical solution, which makes it useful for unit testing.
demes <- ("I")

births <- matrix(c('parms$beta*S*I/(S+I)'), nrow=1, ncol=1)
rownames(births)=colnames(births) <- demes

# migrations undefined for single-deme models
migrations <- matrix(c('0'), nrow=1, ncol=1)
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('parms$mu*I')
names(deaths) <- demes

# non-deme dynamics
ndd <- c(S='(parms$mu)*I - parms$beta*S*I/(S+I)')
