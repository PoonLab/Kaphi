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


update.mstates <- function(z, solve.QAL) {
    z$n.extant <- sum(z$is.extant)

    # process new samples and calculate the state of new lineages
    z$A0 <- z$get.A(z$h0)
    out <- solve.QAL(z$h0, z$h1, z$A0, z$L)
    z$Q <- out[[1]]  # state transition probability matrix
    z$A <- out[[2]]  # update A for this time interval (given Q..)
    z$L <- out[[3]]  # likelihood?

    # clean outputs
    if (is.nan(z$L)) { z$L <- Inf }
    if (any(is.nan(z$Q))) { z$Q <- diag(length(z$A)) }
    if (any(is.nan(z$A))) { z$A <- z$A0 }

    # update mstates (i.e., p_ik(s))
    if (z$n.extant > 1) {
        z$mstates[z$is.extant,] <- t( t(z$Q) %*% t(z$mstates[z$is.extant, ]) )

        z$mstates[z$is.extant, ] <- abs(z$mstates[z$is.extant, ]) /
            rowSums(as.matrix(abs(z$mstates[z$is.extant, ]), nrow=length(z$is.extant)))
        # A_k = sum(p_ik) over i
        z$A <- colSums(as.matrix(z$mstates[z$is.extant, ], nrow=length(z$is.extant)))
    } else {
        # only one extant lineage
        z$mstates[z$is.extant, ] <- t( t(z$Q) %*% z$mstates[z$is.extant, ] )
        z$mstates[z$is.extant, ] <- abs(z$mstates[z$is.extant, ]) / sum(abs(z$mstates[z$is.extant, ]))
        z$A <-  z$mstates[z$is.extant,]  # no need to sum
    }
    return(z)
}



# from R package ECctmc - tweaked to bypass zero row sums check
sample.path <- function (a, b, t0, t1, Q, method = "mr", npaths = 1, eigen_vals = NULL,
    eigen_vecs = NULL, inverse_vecs = NULL, P = NULL)
{
    require(ECctmc)

    if (!method %in% c("mr", "unif")) {
        stop("Simulation method mus be either ", dQuote("mr"),
            " or ", dQuote("unif"))
    }
    if (!all(abs(rowSums(Q) < .Machine$double.eps))) {  # originally "< 0"
        stop("The rate matrix is not valid. The rates must sum to 0 zero within each row.", Q)
    }
    if (!all(diag(Q) <= 0)) {
        stop("The rate matrix is not valid. The diagonal entries must all be non-positive.")
    }
    if (t0 >= t1) {
        stop("t0 must be less than t1.")
    }
    if (!all(c(a, b) %in% 1:nrow(Q))) {
        stop("The endpoints must be given in as row numbers in the rate matrix.")
    }
    if (!is.null(eigen_vals) || !is.null(eigen_vecs) || !is.null(inverse_vecs)) {
        if (is.null(eigen_vals) || is.null(eigen_vecs) || is.null(inverse_vecs)) {
            stop("If one part of the eigen decomposition of Q was provided, all parts must be provided.")
        }
    }
    if (!is.null(P) && !all(rowSums(P) == 1)) {
        stop("A valid transition probability matrix must be provided.")
    }
    if (all(Q[a, ] == 0) & a != b) {
        stop("The process cannot start in an absorbing state if the endpoints are different.")
    }
    if (npaths == 1) {
        if (method == "mr") {
            path <- sample_path_mr(a = a, b = b, t0 = t0, t1 = t1,
                Q = Q)
        }
        else {
            if (is.null(eigen_vals) & is.null(P)) {
                path <- sample_path_unif(a = a, b = b, t0 = t0,
                  t1 = t1, Q = Q)
            }
            else if (!is.null(eigen_vals) & is.null(P)) {
                path <- sample_path_unif2(a = a, b = b, t0 = t0,
                  t1 = t1, Q = Q, eigen_vals = eigen_vals, eigen_vecs = eigen_vecs,
                  inverse_vecs = inverse_vecs)
            }
            else if (is.null(eigen_vals) & !is.null(P)) {
                path <- sample_path_unif3(a = a, b = b, t0 = t0,
                  t1 = t1, Q = Q, P = P)
            }
        }
        colnames(path) <- c("time", "state")
    }
    else {
        path <- vector(mode = "list", length = npaths)
        if (method == "mr") {
            for (k in 1:npaths) {
                path[[k]] <- sample_path_mr(a = a, b = b, t0 = t0,
                  t1 = t1, Q = Q)
                colnames(path[[k]]) <- c("time", "state")
            }
        }
        else {
            if (is.null(P)) {
                P <- comp_expmat(Q = Q * (t1 - t0))
            }
            for (k in 1:npaths) {
                path[[k]] <- sample_path_unif3(a = a, b = b,
                  t0 = t0, t1 = t1, Q = Q, P = P)
                colnames(path[[k]]) <- c("time", "state")
            }
        }
    }
    return(path)
}


coalesce.lineages <- function(z) {
    # current number of sampled lineages at this time point
    z$n.extant <- sum(z$is.extant)

    # retrieve F(s), G(s) and Y(s) for this node height
    this.fgy <- z$get.fgy(z$h1)
    .F <- this.fgy$.F
    .G <- this.fgy$.G
    .Y <- this.fgy$.Y

    a <- z$A / .Y  # normalized lineage counts per deme at coalescent event

    z$extant.lines <- which(z$is.extant)
    .lambdamat <- (t(t(a)) %*% a) * .F  # coalescence hazard

    # pick two demes at random based on number of lineages
    kl <- sample.int(z$m^2, size=1, prob=as.vector(.lambdamat))
    k <- 1 + ((kl-1) %% z$m)  # row
    l <- 1 + floor( (kl-1) / z$m )  # column

    # sample lineages to coalesce from respective demes
    probstates <- as.matrix(z$mstates[z$extant.lines,], nrow=length(z$extant.lines))
    u.i <- sample.int(z$n.extant, size=1, prob=probstates[,k])
    probstates[u.i,] <- 0
    u <- z$extant.lines[u.i]
    v <- sample(z$extant.lines, size=1, prob=probstates[,l])

    z$ustates[u,] <- z$mstates[u,]
    z$ustates[v,] <- z$mstates[v,]
    a.u <- pmin(1, z$mstates[u,] / .Y)
    a.v <- pmin(1, z$mstates[v,] / .Y)

    # coalescence rates by deme states for these two lineages
    lambda.uv <- ( a.u %*% t(a.v) ) * .F + ( a.v %*% t(a.u) ) * .F

    # deme state of new lineage determined by relative proportions of coalescence rates
    palpha <- rowSums(lambda.uv) / sum(lambda.uv)

    # new branch is numbered by lineage counter (lcount)
    alpha <- z$lineage.counter
    z$lineage.counter <- z$lineage.counter + 1

    z$is.extant[alpha] <- TRUE
    z$is.extant[u] <- FALSE  # deactivate lineages that coalesced
    z$is.extant[v] <- FALSE

    z$mstates[alpha,] <- palpha
    z$lstates[alpha,] <- palpha
    z$heights[alpha] <- z$h1

    # update tree variables with new node -- needs alpha, u, v
    z$edge[u,] <- c(alpha, u)
    z$edge.length[u] <- z$h1 - z$heights[u]
    z$edge[v,] <- c(alpha, v)
    z$edge.length[v] <- z$h1 - z$heights[v]


    # construct instantaneous rate matrix for state (deme) transitions
    # see equation (51), Volz Genetics 2012
    if (z$m > 1) {
        qm <- t(a.u * t(.G) + a.u %*% t(1-a/.Y) * t(.F))
        diag(qm) <- -rowSums(qm)
        a.state <- sample(1:z$m, 1, prob=palpha)
        b.state <- sample(1:z$m, 1, prob=z$lstates[u,])
        path <- sample.path(a.state, b.state, 0, z$edge.length[u], qm)
        z$inner.edge[[u]] <- path

        qm <- t(a.v * t(.G) + a.v %*% t(1-a/.Y) * t(.F))
        diag(qm) <- -rowSums(qm)
        b.state <- sample(1:z$m, 1, prob=z$lstates[v,])
        path <- sample.path(a.state, b.state, 0, z$edge.length[v], qm)
        z$inner.edge[[v]] <- path
    }


    return (z)
}


.simulate.ode.tree <- function(sample.times, sample.states, fgy, solve.QAL, A.mx, Et) {
    # Args:
    #   sample.times:
    #   sample.states:
    #   fgy:  List returned from init.fgy()
    #   solve.QAL:  Function returned from init.QAL.solver()
    #   A.mx:  List returned from solve.A.mx()
    #   Et:  List returned from get.event.times()

    z <- list()  # workspace to pass to subfunctions

    # Transfer objects from list arguments
    z$event.times <- Et$event.times
    z$is.sample.event <- Et$is.sample.event
    z$get.fgy <- fgy$func
    z$get.A <- A.mx$get.A

    z$m <- ncol(sample.states)  # number of demes
    z$n <- length(sample.times)  # number of tips (sampled lineages)
    z$S <- 1
    z$L <- 0

    z$max.sample.time <- max(sample.times)
    z$sample.heights <- z$max.sample.time - sample.times

    index <- order(z$sample.heights)
    sorted.sample.heights <- z$sample.heights[index]
    sampled.at.h <- function(h) which(sorted.sample.heights==h)

    sorted.sample.states <- as.matrix(sample.states[index,])

    # initialize variables
    z$Nnode <- z$n-1
    z$num.nodes <- z$Nnode + z$n
    z$edge.length <- rep(-1, z$num.nodes-1)  # does not include root edge
    z$edge <- matrix(-1, nrow=z$num.nodes-1, ncol=2)

    # to cache state transitions along branches
    z$inner.edge <- as.list(1:z$num.nodes-1)


    if (is.null(names(sorted.sample.heights))) {
        z$tip.label <- as.character(1:z$n)  # arbitrary labels
    } else {
        z$tip.label <- names(sorted.sample.heights)
    }

    z$heights <- rep(0, z$num.nodes)
    z$heights[1:z$n] <- sorted.sample.heights

    # matrix of deme states closer to the present
    z$lstates <- matrix(-1, z$num.nodes, z$m)
    z$lstates[1:z$n,] <- sorted.sample.states

    # p_ik in Volz (2012, Genetics):
    #   The probability that branch (i) is in state (k) at time (s) in the past
    z$mstates <- matrix(-1, z$num.nodes, z$m)
    z$mstates[1:z$n,] <- z$lstates[1:z$n,]  # observed sample states

    # matrix of deme states closer to the root
    z$ustates <- matrix(-1, z$num.nodes, z$m)

    # initialize extant lineage statistics with most recent tips (height 0)
    z$h0 <- 0
    z$is.extant <- rep(FALSE, z$num.nodes)
    z$is.extant[sampled.at.h(z$h0)] <- TRUE
    z$extant.lines <- which(z$is.extant)

    # loop over events
    z$lineage.counter <- z$n+1
    for (ih in 1:(length(z$event.times)-1)) {
        z$h0 <- z$event.times[ih]
        z$h1 <- z$event.times[ih+1]

        if (z$is.sample.event[ih+1]) {
            # add new tips
            sat.h1 <- sampled.at.h(h1)
            z$is.extant[sat.h1] <- TRUE
            z$heights[sat.h1] <- z$h1
            next  # continue to next event, bypassing calculations below
        }

		# call helper functions
        z <- update.mstates(z, solve.QAL)
        z <- coalesce.lineages(z)
    }

    # convert tree variables into ape::phylo object
    new.tree <- list(
        edge=z$edge,
        edge.length=z$edge.length,
        Nnode=z$Nnode,
        tip.label=z$tip.label,
        heights=z$heights
    )
    class(new.tree) <- "phylo"
    phylo <- read.tree(text=write.tree(new.tree))

    # reorder edges for compatibility with ape::phylo functions
    sample.times.2 <- sample.times[names(sorted.sample.heights)]
    sample.states.2 <- as.matrix(z$lstates[1:z$n,], nrow=z$n)
    rownames(sample.states.2) <- z$tip.label
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
