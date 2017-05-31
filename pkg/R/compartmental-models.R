require(rcolgem, quietly=TRUE)

##########################################################################################################################
call.rcolgem <- function(nreps, x0, t0, t.end, sampleTimes, sampleStates, births, migrations, deaths, ndd, parms, fgyResolution, integrationMethod) {
  # get numerical solution of ODE system
  fgy <- make.fgy(
          t0,       # start time
          t.end,    # end time
          births,
          deaths,
          ndd,
          x0,
          migrations=migrations,
          parms=parms, 
          fgyResolution=fgyResolution, 
          integrationMethod=integrationMethod
  )
  
  # simulate trees...but which function in rcolgem is it? with '.fgy' or without?
  trees <- simulate.binary.dated.tree.fgy(
    fgy[[1]],  # time axis of ODE solution
    fgy[[2]],  # births
    fgy[[3]],  # migrations
    fgy[[4]],  # deme sizes
    sampleTimes, 
    sampleStates, 
    integrationMethod=integrationMethod, 
    n.reps=nreps   #this isn't a parameter in rcolgem's simulate.binary.dated.tree.fgy
  )
  
  # cast result as a multiPhylo object
  class(trees) <- 'multiPhylo'
  return(trees)
}


#######################################################################################################################
## SIR model w/out vital dynamics, constant population
compartmental.model <- function(theta, nsim, tips, model='si', seed=NA, fgyResolution=500, integrationMethod='adams') {
  '
  Use rcolgem to simulate coalescent trees under susceptible-infected (SI)
  model.
  @param theta : parameter list
  @param nsim : number of replicate trees to simulate
  @param tips : number of tips of zero height (integer) OR vector of tip heights (vector)
  @param model : specified compartmental model ('si', 'sir.nondynamic', 'sir.dynamic', 'sis', 'seir')
  @param seed : set seed for pseudorandom generator
  @param fgyResolution : time resolution of ODE solution
  @param integrationMethod : method for numerical solution of ODE
  '
  # TODO: check contents of theta list
  th.args <- names(theta)
  if (length(th.args) < 4 || any(!is.element(c('t.end', 'N', 'beta', 'gamma', 'mu'), th.args))) {
    stop("'theta' does not hold Kaphi-compatible parameters")
  }

  t0 <- 0  # initial time
  t.end <- theta$t.end  # upper time boundary
  
  
  # initial population frequencies
  S <- theta$N - 1
  I <- 1  # assume epidemic starts with single infected individual
  R <- 0
  x0 <- c(I=I, R=R, S=S)  #sample vector
  
  if (any(x0 < 0)) {
    stop('Population sizes cannot be less than 0.')
  }
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # parsing R expressions representing ODE system
  # can eliminate this and just use theta$ calls
  parms <- list(
    beta=theta$beta,   # transmission rate
    gamma=theta$gamma  # mortality from infection
    mu=theta$mu        # baseline death rate
  )
  if (any(parms < 0)) {
    stop('No negative values permitted for model rate parameters.')
  }
  
  
  ## define ODE system
  
  # demes are subpopulations from which we can sample virus
  demes <- c('I')
  nonDemes <- c('S')
  
  
  # SIR model w/ no vital dynamics
  if (identical(tolower(model), 'sir.nondynamic')) {
    # birth is the rate of lineage splitting - in this case, infection of a susceptible
    # migration is the state transition of a lineage without splitting
    
    # non-deme dynamics is describing the subpopulation
    # note replacement of susceptibles with death of infected, for constant population size
    
    births <- rbind(c('parms$beta * S * I / (S+I)'))
    migrations <- rbind(c('0'))
    deaths <- rbind(c('parms$gamma * I'))
    nonDemeDynamics <- rbind(c('-parms$beta * S * I / (S+I) + parms$gamma * I'))
    
    attr(SIR.nondynamic, 'name') <- 'SIR.nondynamic'  # satisfies requirement in smcConfig.R set.model()
  }
  
  # SIR model w/ births and deaths
  else if (identical(tolower(model), 'sir.dynamic')) {
    births <- rbind(c('parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('0'))
    deaths <- rbind(c('parms$gamma * I - parms$mu * R'))
    nonDemeDynamics <- rbind(c('parms$mu * (I+R) - parms$beta * S * I / (S+I)'))
    
    attr(SIR.dynamic. 'name') <- 'SIR.dynamic'
  }
  
  # SIS model w/ births and deaths
  else if (identical(tolower(model), 'sis')) {
    births <- rbind(c('parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('0'))
    deaths <- rbind(c('parms$gamma * I - parms$mu * R'))
    nonDemeDynamics <- rbind(c('-parms$beta * S * I / (S+I) + parms$mu * (I+R) + parms$gamma * I'))
    
    attr(SIS, 'name') <- 'SIS'
  }
  
  # SEIR model, assuming presence of vital dynamics w/ birth rate equal to death rate
  else if (identical(tolower(model), 'seir')) {
    # first infected will be exposed for an incubation period, not immediately infectious
    # update sample vector to include Exposed compartment
    E <- 1
    I <- 0
    x0 <- c(I=I, R=R, S=S, E=E)
    
    births <- rbind(c('parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('parms$gamma * I - parms$mu * R'))
    nonDemeDynamics <- rbind(c('-parms$beta * S * I / (S+I) + parms$mu * (I+R) + parms$gamma * I'))
    
    # exposed individuals in incubation period
    exposed <- rbind(c('parms$beta * S * I / (S+I) - (parms$episilon + parms$mu) * E'))
    rownames(exposed) <- colnames(exposed) <- demes  # assigned deme state to exposed b/c strictly seir case at the moment
    
    attr(SEIR, 'name') <- 'SEIR'
  }
  
  else {
    stop ("Model is not Kaphi-compatible. Must be a character string of one of the following: 'sir.nondynamic', 'sir.dynamic', 'sis', 'seir' ")
  }
  
  # assigning deme and nondeme states
  rownames(births) <- colnames(births) <- demes
  rownames(migrations) <- colnames(migrations) <- demes
  rownames(deaths) <- colnames(deaths) <- demes 
  names(nonDemeDynamics) <- nonDemes

  
  # sample times
  sampleTimes <- t.end - tip.heights

  # sample states
  sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))
  colnames(sampleStates) <- demes
  rownames(sampleStates) <- 1:n.tips

  
  # calculates numerical solution of ODE system and returns simulated trees
  trees <- call.rcolgem(nsim, x0, t0, t.end, sampleTimes, sampleStates, births, migrations, deaths, ndd, parms, fgyResolution, integrationMethod)
  
  
  return(trees)  # returning an ape phylo object
}  


