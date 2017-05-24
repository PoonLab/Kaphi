## Binary State Speciation and Extinction (BiSSE) Model

require(diversitree, quietly=TRUE)

## Key Probabilities:
##    - DN0(t): Probability that a lineage beginning at time t
##              in state 0 would evolve into a clade like that
##              observed to have descended from node N.
##    - DN1(t): Probability that a lineage beginning at time t
##              in state 1 would evolve into a clade like that
##              observed to have descended from node N.
##    - E0(t): Probability that a lineage beginning at time t
##             in state 0 leaves no descendents at the present day.
##    - E1(t): Probability that a lineage beginning at time t
##             in state 1 leaves no descendents at the present day.

## Parameters:
##    - lambda0: Rate of speciation in state 0.
##    - lambda1: Rate of speciation in state 1.
##    - mu0: Rate of extinction in state 0.
##    - mu1: Rate of extinction in state 1.
##    - q01: Rate of change from state 0 to 1.
##    - q10: Rate of change from state 1 to 0.

#parms <- list(l0=0.1, l1=0.2, m0=0.03, m1=0.03, q01=0.01, q10=0.01)
parms <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)

phy <- tree.bisse(parms, max.taxa=50, max.t=Inf, include.extinct=FALSE, 
                  x0=NA)

states <- phy$tip.state

model <- make.bisse(phy, states)


bisse <- function(theta, nsim, tips, labels=NA, seed=NA) {
    "
    theta : a vector containing the parameter values for the model
    nsim : number of simulations
    tips : either an integer (number of contemporaneous tips) or a
           vector where the length of the vector is number of tips
           and the values correspond to tip heights
              **For BiSSE, ultrametric tree is required**

    Use diversitree to simulate trees under binary state
    speciation and extinction (BiSSE) model.

    @param lambda0 : Rate of speciation in state 0
    @param lambda1 : Rate of speciation in state 1
    @param mu0 : Rate of extinction in state 0
    @param mu1 : Rate of extinction in state 1
    @param q01 : Rate of change from state 0 to 1
    @param q10 : Rate of change from state 1 to 0
    "
  
  if(length(theta) != 6) {
        stop('theta must contain 6 parameters')
  }
  
  if(!is.na(seed)) {
        set.seed(seed)
  }
  
  if(length(tips) < 1) {
        stop('tips must have at least one value')
  } else if(length(tips) > 1) {
        n.tips <- as.integer(length(tips))
  } else {
        n.tips <- as.integer(tips)
  }
  
  #b.trees <- trees(theta, type="bisse", n=nsim, max.taxa=n.tips, max.t=Inf, include.extinct=FALSE)
  tree <- tree.bisse(theta, max.taxa=tips, max.t=Inf, include.extinct=FALSE,
                     x0=NA)
  
  states <- tree$tip.state
  
  model <- make.bisse(tree, states)
  
  #start <- starting.point.bisse(tree)
  
  #fit <- find.mle(lik, p)
  
  #
  
}
