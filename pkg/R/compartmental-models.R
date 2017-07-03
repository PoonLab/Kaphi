# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.


require(rcolgem, quietly=TRUE)


.call.rcolgem <- function(x0, t0, t.end, sampleTimes, sampleStates, births, migrations, deaths, ndd, parms, fgyResolution, integrationMethod) {
  
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
  #plot(fgy[[2]], fgy[[1]])
  #plot(fgy[[3]], fgy[[1]])
  #plot(fgy[[4]], fgy[[1]])
  #plot(fgy[[5]], fgy[[1]])
  #return(fgy)
  
  # simulate tree
  tree <- simulate.binary.dated.tree.fgy(
    fgy[[1]],  # time axis of ODE solution
    fgy[[2]],  # births
    fgy[[3]],  # migrations
    fgy[[4]],  # deme sizes
    sampleTimes,
    sampleStates,
    integrationMethod=integrationMethod  #  , n.reps=nreps   --this isn't a parameter in rcolgem's simulate.binary.dated.tree.fgy
  )
  
  #phy <- tree[[1]]
  #class(phy)  <- 'phylo'
  #return(phy)
  
  # https://github.com/cran/ape/blob/master/R/rtree.R converting tree result into an ape phylo object phy
  phy <- list(edge=tree$edge, edge.length=tree$edge.length)
  if (is.null(tree$tip.label))
    tree$tip.label <- paste("t", 1:tree$n)     # n = number of tips
  phy$tip.label <- sample(tree$tip.label)
  phy$Nnode <- tree$n - 1L
  class(phy) <- "phylo"
  phy <- reorder(phy)
  # to avoid crossings when converting with as.hclust
  phy$edge[phy$edge[,2] <= tree$n, 2] <- 1:tree$n
  phy
  
  return(phy)    # returning an ape phylo object
}


compartmental.model <- function(theta, nsim, tips, model='sir.nondynamic', seed=NA, labels=NA, fgyResolution=500, integrationMethod='adams') {
  "
  Use rcolgem to simulate coalescent trees under susceptible-infected (SI) 
  model.
  @param theta : vector containing parameter values for the model
  @param nsim : number of replicate trees to simulate
  @param tips : number of tips of zero height (integer) OR vector of tip heights (vector)
  @param model : specified compartmental model ('sir.nondynamic', 'sir.dynamic', 'sis', 'seir')
  @param seed : set seed for pseudorandom generator
  @param fgyResolution : time resolution of ODE solution
  @param integrationMethod : method for numerical solution of ODE
  "
  # check contents of theta list
  # params like mu and alpha are only used in certain models
  th.args <- names(theta)
  if (length(th.args) < 6 || any(!is.element(c('t.end', 'N', 'beta', 'gamma', 'mu', 'alpha'), th.args))) {
    stop("'theta' does not hold Kaphi-compatible parameters")
  }

  theta <- as.list(theta)   # convert to list because $ operator is invalid for atomic vectors
  t0 <- 0  # initial time
  t.end <- theta$t.end  # upper time boundary
  
  # initial population frequencies
  S <- theta$N - 1
  I <- 1  # assume epidemic starts with single infected individual
  # x0 <- c(I=I, R=R, S=S)  #sample vector, removed R=R b/c rcolgem allows only births and ndd 1x1 matrices (checks that length(x0) == m + mm)
  x0 <- c(I=I, S=S)
  
  if (any(x0 < 0)) {
    stop('Population sizes cannot be less than 0.')
  }
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # parsing R expressions representing ODE system
  # can eliminate this and just use theta$ calls
  parms <- list(
    beta=theta$beta,        # transmission rate
    gamma=theta$gamma,      # mortality from infection
    mu=theta$mu,            # baseline death rate
    alpha=theta$alpha   # incubation period
  )
  if (any(parms < 0)) {
    stop('No negative values permitted for model rate parameters.')
  }
  
  
  ## define ODE system
  
  # demes are subpopulations from which we can sample virus
  demes <- c('I')
  nonDemes <- c('S')
  
  # SIR model w/ no vital dynamics (births and deaths)
  if (identical(tolower(model), 'sir.nondynamic')) {
    # birth is the rate of lineage splitting - in this case, infection of a susceptible
    # migration is the state transition of a lineage without splitting
    
    # non-deme dynamics is describing the subpopulation
    # note replacement of susceptibles with death of infected, for constant population size
    
    births <- rbind(c('parms$beta * S * I / (S+I) - parms$gamma * I'))
    migrations <- rbind(c('0'))      #in rcolgem manual say this should be omitted if there is only one deme
    deaths <- rbind(c('parms$gamma * I'))
    nonDemeDynamics <- c('-parms$beta * S * I / (S+I)')
  }
  
  # SIR model w/ vital dynamics and constant population
  else if (identical(tolower(model), 'sir.dynamic')) {
    births <- rbind(c('parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('0'))
    deaths <- rbind(c('parms$gamma * I'))          # removed from the population
    nonDemeDynamics <- rbind(c('parms$mu * I - parms$beta * S * I / (S+I)'))
  }
  
  # SIS model w/ births and deaths
  else if (identical(tolower(model), 'sis')) {
    births <- rbind(c('parms$beta * S * I / (S+I) - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('parms$gamma * I'))       # migrating out of I compartment, but back into S compartment
    deaths <- rbind(c('0'))
    nonDemeDynamics <- rbind(c('-parms$beta * S * I / (S+I)'))
  }
  
  # SEIR model, assuming presence of vital dynamics w/ birth rate equal to death rate
  else if (identical(tolower(model), 'seir')) {
    # first infected will be exposed for an incubation period, not immediately infectious
    # update sample vector to include Exposed compartment
    E <- 0
    S <- theta$N - (E+I)
    x0 <- c(S=S, E=E, I=I)
    
    births <- rbind(c('parms$alpha * E - (parms$gamma + parms$mu) * I'))
    migrations <- rbind(c('0'))
    deaths <- rbind(c('parms$gamma * I'))
    nonDemeDynamics <- paste('-parms$beta * S * I / (S+I) + parms$mu * I - parms$mu * E', '+parms$beta * S * I / (S+I) - (parms$alpha + parms$mu) * E', sep='')
    # second part of ndd is the expression for exposed individuals in incubation period
  }
  
  else {
    stop ("Model is not Kaphi-compatible. Must be a character string of one of the following: 'sir.nondynamic', 'sir.dynamic', 'sis', 'seir'.")
  }
  
  # assigning deme and nondeme states
  rownames(births) <- colnames(births) <- demes
  rownames(migrations) <- colnames(migrations) <- demes
  names(deaths) <- demes 
  names(nonDemeDynamics) <- nonDemes

  
  # check if there are tip labels (dates and states) to incorporate
  tips <- .parse.tips(tips)      # function is written in smcConfig.R
  if (sum(x0) < tips$n.tips) {
    stop("Population size is smaller than requested number of tips")
  }
  
  # sample times: a vector indicating when each tip was sampled
  sampleTimes <- sapply(tips$tip.heights, function(x) {t.end - x})

  # sample states
  sampleStates <- matrix(1, nrow=tips$n.tips, ncol=length(demes))
  colnames(sampleStates) <- demes
  rownames(sampleStates) <- 1:tips$n.tips

  
  # calculates numerical solution of ODE system and returns simulated trees
  # incorporate number of simulations
  trees <- replicate(nsim, # num of simulations
                     .call.rcolgem(x0, t0, t.end, sampleTimes, sampleStates, births, migrations=NA, deaths, nonDemeDynamics, parms, fgyResolution, integrationMethod),
                     simplify=FALSE
                     )
  
  # cast result as a multiPhylo object
  class(trees) <- "multiPhylo"
  return(trees)
} 
attr(compartmental.model, 'name') <- "compartmental.model"  # satisfies requirement in smcConfig.R set.model() function

