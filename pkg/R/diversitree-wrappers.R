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
bisse <- function(theta, nsim, tips, labels=NA, seed=NA) {
    "
    theta : a vector containing the parameter values for the model
        @param lambda0 : Rate of speciation in state 0
        @param lambda1 : Rate of speciation in state 1
        @param mu0 : Rate of extinction in state 0
        @param mu1 : Rate of extinction in state 1
        @param q01 : Rate of change from state 0 to 1
        @param q10 : Rate of change from state 1 to 0
    nsim : number of simulations
    tips : if integer, the number of tips; if vector, the height of
           each tip

    Use diversitree to simulate trees under binary state speciation 
    and extinction (BiSSE) model. BiSSE trees are simulated based on
    the rates of speciation and extinction for two character states and
    the rates of changing between these states.
    "
  ## Validate arguments
  if(all(!is.element(c('lambda0', 'lambda1', 'mu0', 'mu1', 'q01', 'q10'), names(theta)))) {
    stop('theta does not contain required parameters')
  }
  ## Set seed
  if(!is.na(seed)) {
        set.seed(seed)
  }
  ## Parse n.tips and tip.heights from tips --> remove/adjust w/ issue #57 resolution
  if(length(tips) < 1) {
        stop('tips must have at least one value')
  } else if(length(tips) > 1) {
        n.tips <- as.integer(length(tips))
        tip.heights <- tips
  } else {
        n.tips <- as.integer(tips)
  }
  ## BiSSE parameter vector
  parms <- unname(theta)
  ## Simulate tree(s)
  result <- trees(parms, type="bisse", n=nsim, max.taxa=n.tips, 
                  max.t=Inf, include.extinct=FALSE)
  return(result)  
}
attr(bisse, 'name') <- "bisse"

## Multiple State Specitation and Extinction (MuSSE) Model
## Quantitative State Speciation and Extinction (QuaSSE) Model
## GeoSSE?
