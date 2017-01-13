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

# calculate at global scope to avoid redundant work
expr <- parse.ode(births, deaths, ndd, migrations)
parms <- list(beta=0.1, mu=0.01)
x0 <- c(S=999, I=1)
sol <- solve.ode(expr, t0=0, t1=100, x0=x0, parms=parms, time.pts=100,
                 integrationMethod='rk4')
