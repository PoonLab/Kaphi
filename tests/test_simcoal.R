require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/SI-matrices.R')

# calculate at global scope to avoid redundant work
expr <- parse.ode(births, deaths, ndd, migrations)
parms <- list(beta=0.1, mu=0.01, gamma=0.02)
x0 <- c(S=999, I=1)
sol <- solve.ode(expr, t0=0, t1=100, x0=x0, parms=parms, time.pts=100,
                 integrationMethod='rk4')


test.init.fgy <- function() {
    result <- init.fgy(sol, max.sample.time=90)
    checkEquals('list', typeof(result))

    fgy.parms <- result$parms
    checkEquals(100, fgy.parms$resolution)
    checkEquals('list', typeof(fgy.parms$F.d))
    checkEquals(100, length(fgy.parms$F.d))
    checkEquals(100, length(fgy.parms$G.d))
    checkEquals(100, length(fgy.parms$Y.d))

    checkTrue(all(fgy.parms$height >= 0))
    checkTrue(all(fgy.parms$height <= 90))
    checkEquals(90, length(fgy.parms$height))

    checkEquals(10, fgy.parms$hoffset)
    checkEquals(1+(100*10)/100, fgy.parms$get.index(0))

    checkEquals(90, nrow(result$mat))  # matrix rows correspond to heights
    checkEquals(1+1+1, ncol(result$mat))  # F, G, and Y are all 1x1 matrices
}

test.solve.A.mx <- function() {
    # note fgy.resolution=100
    fgy <- init.fgy(sol, max.sample.time=90)
    sample.states <- matrix(1, nrow=10, ncol=1)
    sample.heights <- rep(0, times=10)
    result <- solve.A.mx(fgy, sample.states, sample.heights)

    checkEquals(100, nrow(result$mat))
    checkEquals(seq(0, 90, length.out=100), result$haxis)
    checkEquals(10, result$get.A(0))
}


test.init.QAL.solver <- function() {
    fgy <- init.fgy(sol, max.sample.time=90)
    sample.states <- matrix(1, nrow=10, ncol=1)
    sample.heights <- rep(0, times=10)  # all lineages sampled at same time
    result.func <- init.QAL.solver(fgy, sample.states, sample.heights)

    A0 <- 10
    L0 <- 0.1
    # trivial result
    result <- result.func(0, 0, A0, L0)
    checkEquals(as.matrix(1), result[[1]])
    checkEquals(A0, result[[2]])
    checkEquals(L0, result[[3]])

    result <- result.func(0, 1, A0, L0)
    checkTrue(as.vector(result[[1]])>0)
    checkTrue(as.vector(result[[1]])<1)
    checkTrue(result[[2]] < A0)  # number of sampled lineages should decrease back in time
    # TODO: is there a reasonable sanity check on L?
}


test.simulate.ode.tree <- function() {
    sample.times <- rep(90, times=10)  # length determines number of tips
    sample.states <- matrix(1, nrow=10, ncol=1)
    result <- simulate.ode.tree(sol, sample.times, sample.states)
    cat(write.tree(result), "\n")
}
