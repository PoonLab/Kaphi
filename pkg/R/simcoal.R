get.fgymat <- function(sol, t.index, max.sample.time) {
    # ported Erik's code, see if we can refactor
    #  e.g., is it necessary to create an S3 object to compute this matrix?
    times <- sol$times
    fgy.parms <- list(
        resolution=nrow(sol),
        F.d=lapply(t.index, function(k) sol$F[[k]]),
        G.d=lapply(t.index, function(k) sol$G[[k]]),
        Y.d=lapply(t.index, function(k) sol$Y[[k]])
    )
    fgy.parms$hoffset = hoffset <- max(times) - max.sample.time
    if (hoffset < 0) stop("Time axis does not cover the last sample time")
    fgy.parms$get.index <- function(h) {
        min(1 + fgy.parms$resolution * (h + hoffset) / (max(times)-min(times)),
            fgy.parms$resolution)
    }
    fgy.parms$F. <- function(h) { fgy.parms$F.d[[fgy.parms$get.index(h)]] }
    fgy.parms$G. <- function(h) { fgy.parms$G.d[[fgy.parms$get.index(h)]] }
    fgy.parms$Y. <- function(h) { fgy.parms$Y.d[[fgy.parms$get.index(h)]] }
    fgy.parms$h.bounds <- sort(max.sample.time - times)

    get.fgy <- function(h) {
        list(.Y=fgy.parms$Y.(h), .F=fgy.parms$F.(h), .G=fgy.parms$G.(h))
    }

    heights <- sort(fgy.parms$h.bounds)
    fgy.parms$heights <- heights[heights>=0 & heights<=max.sample.time-min(times)]

    fgy.mat <- t(sapply(fgy.parms$heights, function(h)
        with(get.fgy(h), {
            c(as.vector(.F), as.vector(.G), as.vector(.Y))
        })
    ))
    return(pmax(fgy.mat))
}


simulate.ode.tree(sol, sample.times, sample.states, integration.method='rk4') {
    # @param sol:  return value from ode()
    # @param sample.times:  a n-vector of sample collection times, where n is sample size
    # @param sample.states:  an n*m matrix of sample deme states where (m) is number of demes

    ## parse sampleTimes argument
    n.tips <- length(sample.times)
    max.sample.time <- max(sample.times)
    if (is.null(names(sample.times))) {
        # assign arbitrary tip names
        sample.names <- as.character(1:length(sample.times))
        names(sample.times) <- sample.names
        rownames(sample.states) <- sample.names
    }
    sample.heights <- max(sample.times) - sample.times
    ht.index <- order(sample.heights)
    sorted.sample.heights <- sample.heights[ht.index]
    unique.sorted.sample.heights <- unique(sorted.sample.heights)


    ## parse sampleStates argument
    m <- ncol(sample.states)  # number of demes
    deme.names <- colnames(sample.states)
    if (any(!is.element(deme.names, names(sol)))) {
        stop("sampleStates contains deme label not found in ODE solution matrix (sol)")
    }
    if (length(rownames(sample.states)) == 0) {
        stop("sample.states should have row names")
    }
    # cast as matrix to avoid returning a vector in the case of one deme
    sorted.sample.states <- as.matrix(sample.states[ht.index,])


    ## parse times column from ODE solution matrix
    times <- sol$times
    min.time <- min(times)
    max.time <- max(times)
    max.height <- max.sample.time - mintime
    t.index <- order(times, decreasing=TRUE)

    ## construct forcing time series for ODEs
    fgymat <- get.fgymat(sol, t.index, max.sample.time)
}
