get.fgy <- function(sol, max.sample.time) {
    # Args:
    #   sol:  Return value from solve.ode()
    #   t.index:  Rank of ODE time points in descending order
    #   max.sample.time:  max(sample.times)

    # ported Erik's code, see if we can refactor
    #  e.g., is it necessary to create an S3 object to compute this matrix?
    times <- sol$times
    t.index <- order(times, decreasing=TRUE)

    fgy.parms <- list(
        resolution=nrow(sol$sol),
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
    # pmax has no effect with a single argument, I don't know why Erik calls it here...
    return(list(mat=pmax(fgy.mat), func=get.fgy, parms=fgy.parms))
}


init.QAL.solver <- function(fgy, sample.heights) {
    # unpack list argument
    fgy.mat <- fgy$mat
    get.fgy <- fgy$func
    fgy.parms <- fgy$parms

    max.height <- max(sample.heights)

    # internal function to solve for matrices using Erik Volz's C implementation
    function(h0, h1, A0, L0) {
        Q0 <- diag(m)
        parameters <- c(m, max.height, length(fgy.parms$heights), sum(A0), as.vector(fgy.mat))
        y0 <- c(as.vector(Q0), A0, L0)
        o <- ode(y=y0, c(h0, h1), func="dQAL", parms=parameters, dllname="Kaphi",
                 initfunc="initfunc", method=integration.method)
        Q1 <- t(matrix(abs(o[nrow(o), 2:(1+m^2)]), nrow=m))
        A1 <- o[nrow(o), (1+m^2+1):(1+m^2+m)]
        L1 <- o[nrow(o), ncol(o)]
        return(list(unname(Q1), unname(A1), unname(L1)))
    }
}


solve.A.mx <- function(fgy, sample.states, sample.heights) {
    # unpack list argument
    fgy.mat <- fgy$mat
    get.fgy <- fgy$func
    fgy.parms <- fgy$parms

    max.height <- max(sample.heights)
    ht.index <- order(sample.heights)
    sorted.sample.heights <- sample.heights[ht.index]
    unique.sorted.sample.heights <- unique(sorted.sample.heights)

    # cast as matrix to avoid returning a vector in the case of one deme
    sorted.sample.states <- as.matrix(sample.states[ht.index,])

    m <- ncol(sorted.sample.states)  # number of demes
    n <- nrow(sorted.sample.states)  # number of tips

    # cumulative sorted sample/not-sampled states
    cumul.sss <- sapply(1:m, function(k) cumsum(sorted.sample.states[,k]))
    cumul.snss <- t(cumul.sss[n,] - t(cumul.ssss))

    # linear interpolation function on sample heights to locate "not sampled yet" lineages
    nsy.index <- approxfun(
        x=sort(jitter(sorted.sample.heights, factor=max(1e-6, sorted.sample.heights[length(sorted.sample.heights)]/1e6))),
        y=1:n, method='constant', rule=2
    )
    not.sampled.yet <- function(h) {
        cumul.snss[nsy.index(h), ]
    }

    # derivative function for ODE solution
    dA <- function(h, A, parms, ...) {
        nsy <- not.sampled.yet(h)
        with(get.fgy(h), {
            A.Y <- (A-nsy) / .Y
            A.Y[is.nan(A.Y)] <- 0
            csFpG <- colSums(.F+.G)
            list(setNames(as.vector(
                .G %*% A.Y - csFpG * A.Y + (.F %*% A.Y) * pmax(1-A.Y, 0)
            ), names(A)))
        })
    }

    # this appears to be an unused variable!
    A.intervals <- c(sorted.sample.heights[2:length(unique.sorted.sample.heights)], max.height)

    h0 <- 0
    sampled.at.h <- function(h) which(sorted.sample.heights==h)
    haxis <- seq(0, max.height, length.out=fgy.parms$resolution)

    # numerical solution of ODE
    odA <- ode(y=colSums(sorted.sample.states), times=haxis, func=dA, parms=NA, method='adams')

    # unpack results from numerical solution
    haxis <- odA[,1]
    A.plus.unsampled <- odA[,2:(m+1)]
    A.plus.unsampled <- as.matrix(A.plus.unsampled, nrows=length(haxis))

    # return function
    get.A.index <- function(h) {
        min(1 + floor(fgy.parms$resolution * h/max.height), fgy.parms$resolution)
    }
    func <- function(h) {
        i <- get.A.index(h)
        A.plus.unsampled[i,] - not.sampled.yet(h)
    }
    return(list(mat=A.plus.unsampled, haxis=haxis, get.A=func))
}


get.event.times <- function(A.mx, sample.heights) {
    A.plus.unsampled <- A.mx$mat
    haxis <- A.mx$haxis

    A.mono <- rowSums(A.plus.unsampled)
    A.mono[is.na(A.mono)] <- min(A.mono[!is.na(A.mono)])
    A.mono <- A.mono - min(A.mono)
    A.mono <- (max(A.mono)-A.mono) / max(A.mono)

    node.heights <- sort(approx(A.mono, haxis, xout=runif(n-1, 0, 1))$y)
    unique.sample.heights <- unique(sample.heights)
    event.times <- c(unique.sample.heights, node.heights)

    is.sample.event <- c(rep(TRUE, length(unique.sample.heights)), rep(FALSE, length(node.heights)))
    index <- order(event.times)
    return (list(event.times=event.times[index], is.sample.event=is.sample.event[index]))
}


.simulate.ode.tree <- function(sample.times, sample.times) {
    n <- length(sample.times)
    S <- 1
    L <- 0
    max.sample.time <- max(sample.times)
    sample.heights <- max.sample.time - sample.times
    sorted.sample.heights <- sort(sample.heights)

    # initialize variables
    Nnode <- n-1
    edge.length <- rep(-1, Nnode + n-1)  # does not include root edge
    edge <- matrix(-1, nrow=Nnode+n-1, ncol=2)

if (is.null(names(sorted.sample.heights))) {
        tip.label <- as.character(1:n)  # arbitrary labels
    } else {
        tip.label <- names(sorted.sample.heights)
    }

    heights <- rep(0, Nnode+n)
    heights[1:n] <- sorted.sample.heights
    parent.heights <- rep(-1, Nnode+n)
    in.edge.map <- rep(-1, Nnode+n)
    out.edge.map <- matrix(-1, nrow=Nnode+n, ncol=2)
    parent <- 1:(Nnode+n)
    daughters <- matrix(-1, Nnode+n, 2)
}

simulate.ode.tree <- function(sol, sample.times, sample.states, integration.method='rk4') {
    # Args:
    #   sol:  return value from ode()
    #   sample.times:  a n-vector of sample collection times, where n is sample size
    #   sample.states:  an n*m matrix of sample deme states where (m) is number of demes

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

    ## parse sampleStates argument
    m <- ncol(sample.states)  # number of demes
    deme.names <- colnames(sample.states)
    if (any(!is.element(deme.names, colnames(sol$sol)))) {
        stop("sampleStates contains deme label not found in ODE solution matrix (sol)")
    }
    if (length(rownames(sample.states)) == 0) {
        stop("sample.states should have row names")
    }

    ## call helper functions
    fgy <- get.fgy(sol, max.sample.time)
    solve.QAL <- init.QAL.solver(fgy, sample.heights)
    A.mx <- solve.A.mx(fgy, sample.states, sample.heights)
    Et <- get.event.times(A.mx, sample.heights)


}
