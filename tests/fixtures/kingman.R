# pure birth SI model should resemble exponential coalescent
k.births <- matrix('parms$beta*I')
rownames(k.births) <- colnames(births) <- 'I'

k.migrations <- matrix('0')
rownames(k.migrations) <- colnames(k.migrations) <- 'I'

k.deaths <- c(I='0')

k.ndd <- c('-parms$beta*I')
names(k.ndd) <- c('S')

k.expr <- parse.ode(k.births, k.deaths, k.ndd, k.migrations)
k.parms <- list(beta=0.1)
k.x0 <- c(S=999, I=1)
k.sol <- solve.ode(k.expr, t0=0, t1=50, x0=k.x0, parms=k.parms, time.pts=500,
                 integrationMethod='rk4')
