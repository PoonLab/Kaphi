## SIR model w/out vital dynamics, constant population
SIR.nondynamic <- function(theta, nsim, n.tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
  t0 <- 0
  
  # initial population frequencies
  S <- N - 1    # where does the N population parameter come from?
  I <- 1
  R <- 0
  x0 <- c(I=I, S=S, R=R)
  if (any(x0 < 0)) {
    stop("Population sizes cannot be less than 0.")
  }
  
  # parsing R expressions representing ODE system
  parms <- list(beta=beta, gamma=gamma)  # where do the beta and gamma instances come from?
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
  migrations <- rbind(c("parms$gamma * I"))
  rownames(migrations) <- colnames(migrations) <- demes   # demes or a separate one for recovered?
  
  # non-deme dynamics is describing the subpopulation
  nonDemeDynamics <- rbind(c("-parms$beta * S * I / (S+I)"))
  names(nonDemeDynamics) <- nonDemes
  
  return(list(c(births, migrations, nonDemeDynamics)))
}  


######################################################################################################################
## SIR model w/ births and deaths
SIR.dynamic <- function(theta, nsim, n.tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
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
  
  return(list(c(births, migrations, nonDemeDynamics)))
}  

######################################################################################################################
## SIS model with births and deaths
SIS <- function(theta, nsim, n.tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
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
  
  return(list(c(births, nonDemeDynamics)))
}  



######################################################################################################################
## SEIR model, assuming presence of vital dynamics w/ birth rate equal to the death rate
SEIR <- function(theta, nsim, n.tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
  
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
  
  return(list(c(births, migrations, nonDemeDynamics, exposed)))
}  