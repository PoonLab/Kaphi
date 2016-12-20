simulate.ode.tree(sol, sample.times, sample.states, integration.method='rk4') {
    m <- ncol(sampleStates)  # number of demes
    if (any(!is.element(names(sampleStates), names(sol)))) {
        stop("sampleStates contains deme label not found in ODE solution matrix (sol)")
    }

    ## parse sampleTimes argument
    n.tips <- length(sampleTimes)
    if (length(names(sampleTimes)) == 0) {
        # assign arbitrary tip names
        sampleNames <- as.character(1:length(sampleTimes))
        names(sampleTimes) <- sampleNames
        rownames(sampleStates) <- sampleNames
    }
    sampleHeights <- max(sampleTimes) - sampleTimes
    index <- order(sampleHeights)
    sortedSampleHeights <- sampleHeights[index]


    ## parse sampleStates argument
    if (length(rownames(sampleStates)) == 0) {
        stop("sampleStates should have row names")
    }
    # cast as matrix to avoid returning a vector in the case of one deme
    sortedSampleStates <- as.matrix(sampleStates[index,])

}
