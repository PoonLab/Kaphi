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

require(diversitree)

## Binary State Speciation and Extinction (BiSSE) Model
bisse <- function(theta, nsim, n.tips, parms, labels=NA, seed=NA) {
    "
    parms : a vector containing the BiSSE parameters
      @param lambda0 : Rate of speciation in state 0
      @param lambda1 : Rate of speciation in state 1
      @param mu0 : Rate of extinction in state 0
      @param mu1 : Rate of extinction in state 1
      @param q01 : Rate of change from state 0 to 1
      @param q10 : Rate of change from state 1 to 0
    theta : a vector containing the parameter values for the model
    nsim : number of simulations
    tips : either an integer (number of contemporaneous tips) or a
           vector where the length of the vector is number of tips
           and the values correspond to tip heights
              **For BiSSE, ultrametric tree is required**

    Use diversitree to simulate trees under binary state
    speciation and extinction (BiSSE) model.
    "
  if(!is.element('Ne.tau', names(theta))) {
    stop('theta does not contain required parameter "Ne.tau"')
  }
  if(length(parms) != 6) {
   stop('parms requires 6 parameters')
  }
  if(!is.na(seed)) {
        set.seed(seed)
  }
  
  ## For when n.tips is changed to tips (issue #51):
  ## Parse n.tips from tips
  #if(length(tips) < 1) {
  #      stop('tips must have at least one value')
  #} else if(length(tips) > 1) {
  #      n.tips <- as.integer(length(tips))
  #} else {
  #      n.tips <- as.integer(tips)
  #}
  
  result <- lapply(1:3, function(x) {
    tree.bisse(parms, max.taxa=100)
  })
  #result <- lapply(1:nsim, function(x) {
    #tree <- tree.bisse(parms, max.taxa=n.tips, max.t=Inf, 
                       include.extinct=FALSE, x0=NA)
    #tree$edge.length <- tree$edge.length * theta['Ne.tau'] # rescale
    #tree
    #if (!is.na(labels)) {
     # tree$tip.label <- labels
    }
  })
  return(result)  
}
attr(bisse, 'name') <- "bisse"

## Multiple State Specitation and Extinction (MuSSE) Model
## Quantitative State Speciation and Extinction (QuaSSE) Model
## GeoSSE?
