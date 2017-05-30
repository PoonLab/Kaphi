require(rcolgem, quietly=TRUE)


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
  
  # simulate trees
  trees <- simulate.binary.dated.tree(
    fgy[[1]],  # time axis of ODE solution
    fgy[[2]],  # births
    fgy[[3]],  # migrations
    fgy[[4]],  # deme sizes
    sampleTimes, 
    sampleStates, 
    integrationMethod=integrationMethod, 
    n.reps=nreps
  )
  
  # cast result as a multiPhylo object
  class(trees) <- 'multiPhylo'
  return(trees)
}


## SIR model w/out vital dynamics, constant population
SIR.nondynamic <- function(theta, nsim, tips, seed=NA, fgyResolution=500, integrationMethod='adams') {
  "
  Use rcolgem to simulate coalescent trees under susceptible-infected (SI)
  model.
  @param theta : parameter list
  @param nsim : number of replicate trees to simulate
  @param tips : number of tips of zero height (integer) OR vector of tip heights (vector)
  @param seed : set seed for pseudorandom generator
  @param fgyResolution : time resolution of ODE solution
  @param integrationMethod : method for numerical solution of ODE
  "
  # TODO: check contents of theta list

  t0 <- 0  # initial time
  t.end <- theta$t.end
  
  # initial population frequencies
  S <- theta$N - 1

  I <- 1  # assume epidemic starts with single infected individual
  R <- 0
  x0 <- c(I=I, R=R, S=S)
  
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # parsing R expressions representing ODE system
  parms <- list(
    beta=theta$beta,  # transmission rate
    gamma=theta$gamma  # mortality from infection
  )
  if (any(parms < 0)) {
    stop("No negative values permitted for model rate parameters.")
  }
  
  ## define ODE system
  
  # demes are subpopulations from which we can sample virus
  demes <- c("I")
  nonDemes <- c("S")
  
  # birth is the rate of lineage splitting - in this case, infection of a susceptible
  births <- rbind(c("parms$beta * S * I / (S+I)"))
  rownames(births) <- colnames(births) <- demes
  
  # migration is the state transition of a lineage without splitting
  migrations <- rbind(c("0"))
  rownames(migrations) <- colnames(migrations) <- demes
  
  deaths <- rbind(c("parms$gamma * I"))
  rownames(deaths) <- colnames(deaths) <- demes 
  
  # non-deme dynamics is describing the subpopulation
  #  note replacement of susceptibles with death of infected, for constant population size
  nonDemeDynamics <- rbind(c("-parms$beta * S * I / (S+I) + parms$gamma * I"))
  names(nonDemeDynamics) <- nonDemes

  # sample times
  sampleTimes <- t.end - tip.heights

    # sample states
  sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))
  colnames(sampleStates) <- demes
  rownames(sampleStates) <- 1:n.tips
  #maximum t1 (end time)
  max.t.end <- max(sampleTimes)
  
  trees <- call.rcolgem(nsim, x0, t0, t.end, sampleTimes, sampleStates, births, migrations, deaths, ndd, parms, fgyResolution, integrationMethod)
  
  attr(SIR.nondynamic, "name") <- "SIR.nondynamic"  # satisfies requirement in smcConfig.R set.model()
  return(trees)  # returning an ape phylo object
}  


######################################################################################################################
## SIR model w/ births and deaths
SIR.dynamic <- function(theta, nsim, tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {

  t0 <- 0
  
  # initial population frequencies
  S <- N - 1    
  I <- 1
  R <- 0
  x0 <- c(I=I, S=S, R=R)
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  # parsing R expressions representing ODE system
  parms <- list(beta=beta, gamma=gamma, mu=mu)  
  if (any(parms < 0)) {
    stop("No negative values permitted for model rate parameters.")
  }
  
  # define ODE system
  
  demes <- c("I")
  nonDemes <- c("S")
  
  births <- rbind(c("parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I"))
  rownames(births) <- colnames(births) <- demes
  
  migrations <- rbind(c("parms$gamma * I - parms$mu * R"))
  rownames(migrations) <- colnames(migrations) <- demes 

  nonDemeDynamics <- rbind(c("parms$mu * (I+R) - parms$beta * S * I / (S+I)"))
  names(nonDemeDynamics) <- nonDemes
  
  attr(SIR.dynamic. "name") <- "SIR.dynamic"
  
  return(list(c(births, migrations, nonDemeDynamics)))
}  

######################################################################################################################
## SIS model with births and deaths
SIS <- function(theta, nsim, tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {

  t0 <- 0
  
  # initial population frequencies
  S <- N - 1    
  I <- 1
  R <- 0
  x0 <- c(I=I, S=S, R=R)
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  # parsing R expressions representing ODE system
  parms <- list(beta=beta, gamma=gamma, mu=mu)  
  if (any(parms < 0)) {
    stop("No negative values permitted for model rate parameters.")
  }
  
  # define ODE system
  
  demes <- c("I")
  nonDemes <- c("S")
  
  births <- rbind(c("parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I"))
  rownames(births) <- colnames(births) <- demes
  
  nonDemeDynamics <- rbind(c("-parms$beta * S * I / (S+I) + parms$mu * (I+R) + parms$gamma * I"))
  names(nonDemeDynamics) <- nonDemes
  
  attr(SIS, "name") <- "SIS"
  
  return(list(c(births, nonDemeDynamics)))
}  



######################################################################################################################
## SEIR model, assuming presence of vital dynamics w/ birth rate equal to the death rate
SEIR <- function(theta, nsim, tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {

  t0 <- 0
  
  # initial population frequencies
  S <- N - 1    
  I <- 1
  R <- 0
  E <- 0
  x0 <- c(I=I, S=S, R=R, E=E)
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  # parsing R expressions representing ODE system
  parms <- list(beta=beta, gamma=gamma, mu=mu, epsilon=epsilon)   # epsilon is the incubation period 
  if (any(parms < 0)) {
    stop("No negative values permitted for model rate parameters.")
  }
  
  # define ODE system
  
  demes <- c("I")
  nonDemes <- c("S")
  exp <- c("E")
  
  births <- rbind(c("parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I"))
  rownames(births) <- colnames(births) <- demes
  
  migrations <- rbind(c("parms$gamma * I - parms$mu * R"))
  rownames(migrations) <- colnames(migrations) <- demes
  
  nonDemeDynamics <- rbind(c("-parms$beta * S * I / (S+I) + parms$mu * (I+R) + parms$gamma * I"))
  names(nonDemeDynamics) <- nonDemes
  
  # exposed individuals in incubation period
  exposed <- rbind(c("parms$beta * S * I / (S+I) - (parms$episilon + parms$mu) * E"))
  rownames(exposed) <- colnames(exposed) <- exp
  
  attr(SEIR, "name") <- "SEIR"
  
  return(list(c(births, migrations, nonDemeDynamics, exposed)))
}
