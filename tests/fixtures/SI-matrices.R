# rcolgem definition for a standard SI model
demes <- ("I")

births <- matrix(c('parms$beta*S*I/(S+I)'), nrow=1, ncol=1)
rownames(births)=colnames(births) <- demes

migrations <- matrix(0, nrow=1, ncol=1)
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('(parms$mu+parms$gamma)*I')
names(deaths) <- demes

# non-deme dynamics
ndd <- c(S='(parms$mu+parms$gamma)*I - parms$beta*S*I/(S+I)')
