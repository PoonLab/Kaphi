require(rcolgem, quietly=TRUE)


## SIR model w/out vital dynamics, constant population
SIR.nondynamic <- function(theta, nsim, tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
  }
  
  "
  rcolgem is used to simulate coalescent trees under susceptible-infected-recovered (SIR) model
  @param t.end : boundary condition for ODE solution, time scale of simulation
  @param N : total population size (constant)
  @param beta : transmission rate 
  @param gamma : additional mortality from infection
  @param mu : baseline mortality rate
  @param fgyResolution : time resolution of ODE solution
  @param integrationMethod : method for numerical solution of ODE
  "
  
  t0 <- 0
  t.end <- 30.*52 # need some kind of end time initialized
  
  # initial population frequencies
  S <- N - 1    # where does the N population parameter come from? Needs to be initialized or passed through formals
  I <- 1
  R <- 0
  x0 <- c(I=I, R=R, S=S)
  
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # parsing R expressions representing ODE system
  parms <- list(beta=beta, gamma=gamma)  # where do the beta and gamma instances come from? # initializes model parameters in kamphir/drivers/simulate.SI.R
  if (any(parms < 0)) {
    stop("No negative values permitted for model rate parameters.")
  }
  
  # define ODE system
  
  # demes are subpopulations from which we can sample virus
  demes <- c("I")
  nonDemes <- c("S")
  
  #birth is the rate of lineage splitting - in this case, infection of a susceptible
  births <- rbind(c("parms$beta * S * I / (S+I) - parms$gamma * I"))
  rownames(births) <- colnames(births) <- demes
  
  # migration is the state transition of a lineage without splitting
  migrations <- rbind(c("0"))
  rownames(migrations) <- colnames(migrations) <- demes
  
  deaths <- rbind(c("parms$gamma * I"))
  rownames(deaths) <- colnames(deaths) <- demes 
  
  # non-deme dynamics is describing the subpopulation
  nonDemeDynamics <- rbind(c("-parms$beta * S * I / (S+I)"))
  names(nonDemeDynamics) <- nonDemes
  
  
  # sample times
  sampleTimes <- t.end - tip.heights   #you want tip heights extracted from somewhere, or initialized
  # sample states
  sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))
  colnames(sampleStates) <- demes
  rownames(sampleStates) <- 1:n.tips
  #maximum t1 (end time)
  max.t.end <- max(sampleTimes)
  
  
  #numerical solution of ODE, calling rcolgem function
  tode <- make.fgy(t0, max.t.end, births, deaths, nonDemeDynamics, x0, migrations=migrations, parms=parms, time.pts=2000, integrationMethod=integrationMethod)
  
  #simulating trees, calling rcolgem function
  trees <- simulate.binary.dated.tree.fgy(tode[[2]], tode[[3]], tode[[4]], tode[[1]], x0, sampleTimes, sampleStates, migrations=migrations, parms=parms, integrationMethod=integrationMethod, n.reps=nsim)
  
  attr(SIR.nondynamic, "name") <- "SIR.nondynamic" # satisfies requirement in smcConfig.R set.model()
  return(trees) # returning an ape phylo object
}  


######################################################################################################################
## SIR model w/ births and deaths
SIR.dynamic <- function(theta, nsim, tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
  }
  
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
  
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
  }
  
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
  
  if(length(tips) < 1) {
    stop('tips must have at least one value')
  } else if(length(tips) > 1) {
    n.tips <- as.integer(length(tips))
    tip.heights <- tips
  } else {
    n.tips <- as.integer(tips)
  }
  
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