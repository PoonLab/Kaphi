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

init.fgy <- function(sol, max.sample.time) {
    # Rescale the time axis of F/G/Y matrices in reverse time with maximum
    # sample time as origin.
    #
    # Args:
    #   sol:  Return value from solve.ode()
    #   max.sample.time:  max(sample.times)
    #
    # Returns:
    #   A list containing:
    #     mat: a matrix of F, G, and Y matrix rows that are concatenated into a
    #          single vector per time point;
    #     func: a reference to the function used to generate the matrices
    #     parms: the list of parameters used in calculating the matrices
    times <- sol$times
    t.index <- order(times, decreasing=TRUE)

    fgy.parms <- list(
        resolution=nrow(sol$sol),
        F.d=lapply(t.index, function(k) sol$F[[k]]),
        G.d=lapply(t.index, function(k) sol$G[[k]]),
        Y.d=lapply(t.index, function(k) sol$Y[[k]])
    )

    # reverse time scale with max.sample.time as origin
    heights <- sort(max.sample.time - times)
    fgy.parms$heights <- heights[heights>=0 & heights<=max.sample.time-min(times)]

    # time difference between last sample time and right bound of ODE solution (t1)
    fgy.parms$hoffset = hoffset <- max(times) - max.sample.time
    if (hoffset < 0) stop("Time axis does not cover the last sample time")

    # convert height argument to index on time-axis of ODE solution
    #  at h = max.sample.time-min(times), returns max index (fgy.parms$resolution)
    #  at h = 0, returns
    fgy.parms$get.index <- function(h) {
        min(1 + fgy.parms$resolution * (h + hoffset) / (max(times)-min(times)),
            fgy.parms$resolution)
    }
    fgy.parms$F. <- function(h) { fgy.parms$F.d[[fgy.parms$get.index(h)]] }
    fgy.parms$G. <- function(h) { fgy.parms$G.d[[fgy.parms$get.index(h)]] }
    fgy.parms$Y. <- function(h) { fgy.parms$Y.d[[fgy.parms$get.index(h)]] }

    get.fgy <- function(h) {
        list(.Y=fgy.parms$Y.(h), .F=fgy.parms$F.(h), .G=fgy.parms$G.(h))
    }

    fgy.mat <- t(sapply(fgy.parms$heights, function(h)
        with(get.fgy(h), {
            c(as.vector(.F), as.vector(.G), as.vector(.Y))
        })
    ))
    return(list(mat=fgy.mat, func=get.fgy, parms=fgy.parms))
}


init.QAL.solver <- function(fgy, sample.states, sample.heights, integration.method='rk4') {
    # Q is the state transition (deme migration) rate matrix
    # A is the number of extant (sampled) lineages over time
    # L is the cumulative hazard of coalescence
    #
    # Args:
    #   fgy:  List returned from get.fgy()
    #   sample.heights:  A vector of sampling times relative to most recent
    #                    sample.
    # Returns:
    #   Reference to an internal function for computing Q, A, and L1 matrices
    fgy.mat <- fgy$mat
    get.fgy <- fgy$func
    fgy.parms <- fgy$parms

    # passed along in returned function's environment
    m <- ncol(sample.states)  # number of demes
    heights <- fgy.parms$heights
    max.height <- max(heights)

    # internal function to solve for matrices using Erik Volz's C implementation
    function(h0, h1, A0, L0) {
        Q0 <- diag(m)
        parameters <- c(m, max.height, length(fgy.parms$heights), sum(A0), as.vector(fgy.mat))

        # y is a concatenation of vectorized matrices Q, A, and L
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
    # The A matrix represents the expected number of extant (sampled) lineages over time.
    # Args:
    #   fgy:  List returned from get.fgy()
    #   sample.states:
    #   sample.heights:
    # Returns:
    #   A list containing numerical solution for A matrix, time axis and an accessor function

    # unpack list argument
    fgy.mat <- fgy$mat
    get.fgy <- fgy$func
    fgy.parms <- fgy$parms

    heights <- fgy.parms$heights
    max.height <- max(heights)  # max.sample.time - min(times)

    ht.index <- order(sample.heights)
    sorted.sample.heights <- sample.heights[ht.index]
    unique.sorted.sample.heights <- unique(sorted.sample.heights)

    # cast as matrix to avoid returning a vector in the case of one deme
    sorted.sample.states <- as.matrix(sample.states[ht.index,])

    m <- ncol(sorted.sample.states)  # number of demes

    not.sampled.yet <- function(h) {
        uniq.index <- as.integer(cut(h, breaks=c(unique.sorted.sample.heights, Inf), right=FALSE))
        return(sum(sorted.sample.heights > unique.sorted.sample.heights[uniq.index]))
    }

    # derivative function for ODE solution - see equation (56) from Volz (2012, Genetics)
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

    # numerical solution of ODE
    haxis <- seq(0, max.height, length.out=fgy.parms$resolution)
    odA <- ode(y=colSums(sorted.sample.states), times=haxis, func=dA, parms=NA, method='adams')

    # unpack results from numerical solution
    haxis <- odA[,1]
    A.plus.unsampled <- odA[,2:(m+1)]  # both sampled and not-yet-sampled lineages over time
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
    #  An event is either:
    #   1. sampling of a lineage
    #   2. coalescence of two lineages or
    #   3. a state transition where a single lineage migrates between demes.
    #  Currently this function only looks at the first two event types.

    A.plus.unsampled <- A.mx$mat
    haxis <- A.mx$haxis
    n <- length(sample.heights)  # number of tips

    A.mono <- rowSums(A.plus.unsampled)  # sum across demes per time point
    A.mono[is.na(A.mono)] <- min(A.mono[!is.na(A.mono)])  # impute missing entries
    A.mono <- A.mono - min(A.mono)  # shift values to zero lower bound
    A.mono <- (max(A.mono)-A.mono) / max(A.mono)  # invert and normalize to range [0,1]

    # note max(A.mono) should coincide with most recent time point

    # sample n-1 random points along normalized A.mono trajectory
    # use interpolation to impute heights for these internal nodes
    node.heights <- sort(approx(x=A.mono, y=haxis, xout=runif(n-1, 0, 1))$y)
    unique.sample.heights <- unique(sample.heights)

    # events occur at internal nodes of the tree
    event.times <- c(unique.sample.heights, node.heights)

    is.sample.event <- c(rep(TRUE, length(unique.sample.heights)), rep(FALSE, length(node.heights)))
    index <- order(event.times)
    return (list(event.times=event.times[index], is.sample.event=is.sample.event[index]))
}


.simulate.ode.tree <- function(sample.times, sample.states, fgy, solve.QAL, A.mx, Et) {
    # Args:
    #   sample.times:
    #   sample.states:
    #   fgy:  List returned from init.fgy()
    #   solve.QAL:  Function returned from init.QAL.solver()
    #   A.mx:  List returned from solve.A.mx()
    #   Et:  List returned from get.event.times()

    # Unpack objects from list arguments
    event.times <- Et$event.times
    is.sample.event <- Et$is.sample.event
    get.fgy <- fgy$func
    get.A <- A.mx$get.A

    m <- ncol(sample.states)  # number of demes
    n <- length(sample.times)  # number of tips (sampled lineages)
    S <- 1
    L <- 0

    max.sample.time <- max(sample.times)
    sample.heights <- max.sample.time - sample.times
    index <- order(sample.heights)
    sorted.sample.heights <- sample.heights[index]
    sampled.at.h <- function(h) which(sorted.sample.heights==h)

    sorted.sample.states <- as.matrix(sample.states[index,])

    # initialize variables
    Nnode <- n-1
    num.nodes <- Nnode + n
    edge.length <- rep(-1, Nnode + n-1)  # does not include root edge
    edge <- matrix(-1, nrow=num.nodes-1, ncol=2)

    if (is.null(names(sorted.sample.heights))) {
        tip.label <- as.character(1:n)  # arbitrary labels
    } else {
        tip.label <- names(sorted.sample.heights)
    }

    heights <- rep(0, num.nodes)
    heights[1:n] <- sorted.sample.heights

    # matrix of deme states closer to the present
    lstates <- matrix(-1, Nnode+n, m)
    lstates[1:n,] <- sorted.sample.states

    # p_ik in Volz (2012, Genetics):
    #   The probability that branch (i) is in state (k) at time (s) in the past
    mstates <- matrix(-1, Nnode+n, m)
    mstates[1:n,] <- lstates[1:n]  # observed sample states

    # matrix of deme states closer to the root
    ustates <- matrix(-1, Nnode+n, m)

    # initialize extant lineage statistics with most recent tips (height 0)
    h0 <- 0
    is.extant <- rep(FALSE, Nnode+n)
    is.extant[sampled.at.h(h0)] <- TRUE
    extant.lines <- which(is.extant)

    # A0 is a vector of the number of lineages at height 0 per deme
    if (length(extant.lines) > 1) {
        A0 <- colSums(as.matrix(sorted.sample.states[extant.lines, ], nrow=length(extant.lines)))
    } else {
        # handle the case of a single deme
        A0 <- sorted.sample.states[extant.lines, ]
    }

    # loop over events
    lineage.counter <- n+1
    for (ih in 1:(length(event.times)-1)) {
        h0 <- event.times[ih]
        h1 <- event.times[ih+1]
        fgy <- get.fgy(h1)

        # current number of sampled lineages at this time point
        n.extant <- sum(is.extant)

        # process new samples and calculate the state of new lineages
        A0 <- get.A(h0)
        out <- solve.QAL(h0, h1, A0, L)
        Q <- out[[1]]  # state transition probability matrix
        A <- out[[2]]  # update A for this time interval (given Q..)
        L <- out[[3]]  # likelihood?  It's not used here...

        # clean outputs
        if (is.nan(L)) { L <- Inf }
        if (any(is.nan(Q))) { Q <- diag(length(A)) }
        if (any(is.nan(A))) { A <- A0 }

        # update mstates (i.e., p_ik(s))
        if (n.extant > 1) {
            mstates[is.extant, ] <- t( t(Q) %*% mstates[is.extant, ] )
            mstates[is.extant, ] <- abs(mstates[is.extant, ]) /
                rowSums(as.matrix(abs(mstates[is.extant, ]), nrow=length(is.extant)))
            # A_k = sum(p_ik) over i
            A <- colSums(as.matrix(mstates[is.extant, ], nrow=length(is.extant)))
        } else {
            mstates[is.extant, ] <- t( t(Q) %*% mstates[is.extant, ] )
            mstates[is.extant, ] <- abs(mstates[is.extant, ]) / sum(abs(mstates[is.extant, ]))
            A <-  mstates[is.extant,]  # only one extant lineage, no summation necessary
        }

        if (is.sample.event[ih+1]) {
            # add new tips
            sat.h1 <- sampled.at.h(h1)
            is.extant[sat.h1] <- TRUE
            heights[sat.h1] <- h1
        } else {


            # coalescent event
            .F <- fgy$.F
            .G <- fgy$.G
            .Y <- fgy$.Y

            a <- A / .Y  # normalized lineages over time

            extant.lines <- which(is.extant)
            .lambdamat <- (t(t(a)) %*% a) * .F  # coalescence hazard

            # pick two demes at random based on number of lineages
            kl <- sample.int(m^2, size=1, prob=as.vector(.lambdamat))
            k <- 1 + ((kl-1) %% m)  # row
            l <- 1 + floor( (kl-1) / m )  # column

            # sample lineages to coalesce from respective demes
            probstates <- as.matrix(mstates[extant.lines,], nrow=length(extant.lines))
            u.i <- sample.int(n.extant, size=1, prob=probstates[,k])
            probstates[u.i,] <- 0
            u <- extant.lines[u.i]
            v <- sample(extant.lines, size=1, prob=probstates[,l])

            ustates[u,] <- mstates[u,]
            ustates[v,] <- mstates[v,]
            a.u <- pmin(1, mstates[u,] / .Y)
            a.v <- pmin(1, mstates[v,] / .Y)

            # matrix of coalescence rates by deme states
            lambda.uv <- ( a.u %*% t(a.v) ) * .F + ( a.v %*% t(a.u) ) * .F

            # deme state of new lineage determined by relative proportions of coalescence rates
            palpha <- rowSums(lambda.uv) / sum(lambda.uv)

            # new branch is numbered by counter
            alpha <- lineage.counter
            lineage.counter <- lineage.counter + 1
            is.extant[alpha] <- TRUE
            is.extant[u] <- FALSE  # deactivate lineages that coalesced
            is.extant[v] <- FALSE

            mstates[alpha,] <- palpha
            lstates[alpha,] = palpha  # why use this assignment operator?
            heights[alpha] <- h1

            # update tree variables with new node
            edge[u,] <- c(alpha, u)
            edge.length[u] <- h1 - heights[u]
            edge[v,] <- c(alpha, v)
            edge.length[v] <- h1 - heights[v]
        }
    }

    # convert tree variables into ape::phylo object
    new.tree <- list(
        edge=edge,
        edge.length=edge.length,
        Nnode=Nnode,
        tip.label=tip.label,
        heights=heights
    )
    class(new.tree) <- "phylo"
    phylo <- read.tree(text=write.tree(new.tree))

    # reorder edges for compatibility with ape::phylo functions
    sample.times.2 <- sample.times[names(sorted.sample.heights)]
    sample.states.2 <- as.matrix(lstates[1:n,], nrow=n)
    rownames(sample.states.2) <- tip.label
    sample.times.2 <- sample.times.2[phylo$tip.label]
    sample.states.2 <- as.matrix(sample.states.2[phylo$tip.label, ], nrow=length(phylo$tip.label))

    phylo$sampleTimes <- sample.times.2
    phylo$sampleStates <- sample.states.2

    return(phylo)
}



simulate.ode.tree <- function(sol, sample.times, sample.states, integration.method='rk4') {
    # Args:
    #   sol:  return value from ode()
    #   sample.times:  a n-vector of sample collection times, where n is sample size
    #   sample.states:  an n*m matrix of sample deme states where (m) is number of demes
    #
    # Returns:
    #   A tree (ape phylo object)

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
    fgy <- init.fgy(sol, max.sample.time)
    solve.QAL <- init.QAL.solver(fgy, sample.states, sample.heights)
    A.mx <- solve.A.mx(fgy, sample.states, sample.heights)
    Et <- get.event.times(A.mx, sample.heights)

    sim.tree <- .simulate.ode.tree(sample.times, sample.states, fgy, solve.QAL, A.mx, Et)
    return(sim.tree)
}
