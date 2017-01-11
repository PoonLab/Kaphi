require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/SI-matrices.R')

test.init.fgy <- function() {
    expr <- parse.ode(births, deaths, ndd, migrations)
    parms <- list(beta=0.1, mu=0.01, gamma=0.02)
    x0 <- c(S=999, I=1)
    sol <- solve.ode(expr, t0=0, t1=100, x0=x0, parms=parms, time.pts=100,
                     integrationMethod='rk4')
    result <- init.fgy(sol, max.sample.time=90)

    checkEquals('list', typeof(result))
    checkEquals(100, nrow(result$fgy.mat))
    checkEquals(3, ncol(result$fgy.mat))
}