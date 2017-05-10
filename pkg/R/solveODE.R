require(deSolve)

# Based on Erik Volz's implementation in rcolgem (https://github.com/rforge/colgem)

# @article{volz2012complex,
#   title={Complex population dynamics and the coalescent under neutrality},
#   author={Volz, Erik M},
#   journal={Genetics},
#   volume={190},
#   number={1},
#   pages={187--201},
#   year={2012},
#   publisher={Genetics Soc America}
# }


# Demes are subpopulations from which lineages can be sampled.

# Non-demes are subpopulations from which lineages are NOT sampled, but
#   contribute to the population growth process.

parse.ode <- function(births, deaths, ndd, migrations) {
    # parse string matrices into unevaluated R expressions
    result <- list(
        demeNames=rownames(births),
        nonDemeNames=names(ndd),
        pbirths = parse(text=t(births)),  # transpose to return row-wise
        pmigrations = parse(text=t(migrations)),
        pdeaths = parse(text=deaths),
        pndd = ifelse(!is.na(ndd) & length(ndd)>0, parse(text=ndd), NA)
    )
    return(result)
}


.get.times <- function(t0, t1, time.pts) {
    # generate vector of time points
    times1 <- seq(t0, t1, length.out=time.pts)
    tdelta <- times1[2] - times1[1]  # time difference

    # not clear to me why this is necessary...
    if ( (t0-tdelta)<=t0 ) {
        times0 <- c()
    } else {
        times0 <- seq(t0, t0-tdelta, by=tdelta)
    }
    times <- unique(c(times0, times1))
    return(times)
}



solve.ode <- function(expr, t0, t1, x0, parms, time.pts=2000, integrationMethod='rk4') {
    # Numerical solution of ODE system
    #  derived from rcolgem:make.fgy()
    #
    # Args:
    #   expr:  Parsed expressions from parse.ode()
    #   t0:  Initial time for ODE solution
    #   t1:  End time for ODE solution
    #   x0:  Initial state vector
    #   parms:  Model parameter settings
    #   time.pts:  Number of time points in range [t0, t1]
    #   integrationMethod:  passed to ode()
    #
    # Returns a list containing (K = time.pts):
    #   times:  (K) time points of numerical solutions in reverse order
    #   Y:  K (m) vectors of deme population sizes
    #   F:  K (m x m) matrices of birth rates
    #   G:  K (m x m) matrices of migration rates
    #   sol:  Return value of ode()

    # unpack expressions list
    demeNames <- expr$demeNames
    nonDemeNames <- expr$nonDemeNames
    m <- length(demeNames)      # number of demes
    mm <- length(nonDemeNames)  # number of non-demes

    # validate initial conditions argument
    if (length(x0) != m+mm) {
        stop("Initial conditions incorrect dimension", x0, m, mm)
    }
    if ( any(!is.element(c(demeNames, nonDemeNames), names(x0))) ) {
        stop("Initial conditions vector is incomplete", names(x0), demeNames, nonDemeNames)
    }
    x0 <- x0[c(demeNames, nonDemeNames)]  # filter and sort initial conditions vector

    # process time arguments
    if (time.pts < 1) stop("time.pts argument must be >= 1")
    if (t0 > t1) stop("Start time t0 must be less than end time t1")
    times <- .get.times(t0, t1, time.pts)


    # internal functions
    .birth.matrix <- function(x, t) {
        with(as.list(x),
            t(matrix(sapply(1:m^2, function(k) eval(expr$pbirths[k])),
                     nrow=m, ncol=m))
        )
    }
    .migration.matrix <- function(x, t) {
        with(as.list(x),
            t(matrix(sapply(1:m^2, function(k) eval(expr$pmigrations[k])),
                     nrow=m, ncol=m))
        )
    }
    tBirths <- function(x, t) { colSums(.birth.matrix(x,t)) }
    tMigrationsIn <- function(x, t) { colSums(.migration.matrix(x,t)) }
    tMigrationsOut <- function(x, t) { rowSums(.migration.matrix(x,t)) }
    tDeaths <- function(x, t) {
        with(as.list(x, t), sapply(1:m, function(k) eval(expr$pdeaths[k])))
    }
    dNonDeme <- function(x, t) {
        with(as.list(x, t), sapply(1:mm, function(k) eval(expr$pndd[k])))
    }


    # Function that computes values of derivatives in ODE system (see help(ode))
    dx <- function(t, y, parms, ...) {
        # see equation (28) in Volz 2012 Genetics
        dxdeme <- setNames(
            tBirths(y,t) + tMigrationsIn(y,t) - tMigrationsOut(y,t) - tDeaths(y,t),
        demeNames)
        dxnondeme <- ifelse(mm>0, setNames(dNonDeme(y,t), nonDemeNames), NULL)
        list(c(dxdeme, dxnondeme))  # return values
    }

    # numerical solution of ODE
    sol <- ode(y=x0, times=times, func=dx, parms=parms, method=integrationMethod)

    # prepare return values
    index <- nrow(sol):1
    result <- list(
        times=rev(times),
        Y=lapply(index, function(i) sol[i, demeNames]),
        F=lapply(index, function(i) .birth.matrix(sol[i,], sol[i,1])),
        G=lapply(index, function(i) .migration.matrix(sol[i,], sol[i,1])),
        sol=sol
    )
    return(result)
}
