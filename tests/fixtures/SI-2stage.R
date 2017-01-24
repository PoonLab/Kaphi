# rcolgem definition for an SI model with acute and chronic stages of infection
demes <- c("I1", "I2")

births <- matrix(
  c('parms$beta1*S*I1 / (S+I1+I2)', '0',   # S + I1 -> I1 + I1
    'parms$beta2*S*I2 / (S+I1+I2)', '0'),  # S + I2 -> I1 + I2
nrow=2, byrow=T)
rownames(births)<-colnames(births) <- demes

# migrations undefined for single-deme models
migrations <- matrix(
    c('0', 'parms$alpha * I1',  # I1 -> I2
      '0', '0'),
nrow=2, byrow=TRUE)
rownames(migrations)<-colnames(migrations) <- demes

deaths <- c('(parms$mu)*I1', '(parms$mu+parms$gamma)*I2')
names(deaths) <- demes

# non-deme dynamics
ndd <- c(S='parms$mu*I1 + (parms$mu+parms$gamma)*I2 - S*(parms$beta1*I1 + parms$beta2*I2) / (S+I1+I2)')

# calculate at global scope to avoid redundant work
expr <- parse.ode(births, deaths, ndd, migrations)
parms <- list(beta1=0.01, beta2=0.001, alpha=0.01, gamma=1/520., mu=1/3640.)
x0 <- c(S=999, I1=1, I2=0)
sol <- solve.ode(expr, t0=0, t1=30.*52, x0=x0, parms=parms, time.pts=500,
                 integrationMethod='rk4')
