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

## Standard Speciation Models
speciation.model <- function(theta, nsim, tips, model, seed=NA, ...) {
    "
    theta : a vector containing the parameter values for the model
        @param lambda(i) : Rate of speciation (for state i)
        @param lambda(ijk) : Rate of speciation from ancestral state i to
                             daughter states j and k.
        @param mu(i) : Rate of extinction (for state i)
        @param qij : Rate of change from state i to j
        @param pic : Probabilitiy that a character change from ancetral
                     state i occures during speciation
        @param pia : Probability that a character change from ancestral
                     state i during speciation is asymetrical
        @param sA : Rate of speciation within region A
        @param xA : Rate of extinction/extirpation from region A.
        @param dA : Rate of dispersal/range expansion from region A.
        @param char : model of character evolution.
    nsim : number of trees to simulate
    tips : if integer, the number of tips; if vector, the heights of the tips
    model : model under which trees are evolved. One of:
              - bisse       -- Binary State Speciation and Extinction
              - bisseness   -- BiSSE-Node Enhanced State Shift
              - bd          -- Birth-Death model
              - classe      -- Cladogenetic State change Speciation and Extinction
              - geosse      -- Geographic State Speciation and Extinction
              - musse       -- Multi-state Speciation and Extinction
              - quasse      -- Quantitative state Specatiation and Extinction
              - yule        -- Yule model

    Use diversitree to simulate trees under one of the 'SSE' models or a simple
    character indepdenent Birth-Death or Yule model.
    "
  ## Validate arguments
  if(model == 'bisse') {
    if(any(!is.element(c('lambda0', 'lambda1', 'mu0', 'mu1', 'q01', 'q10'),
                       names(theta)))) {
      stop('theta does not contain all of the required parameters:
           lambda0, lambda1, mu0, mu1, q01, q10')
    }
  } else if(model == 'bisseness') {
      if(any(!is.element(c('lambda0', 'lambda1', 'mu0', 'mu1', 'q01', 'q10',
                           'p0c', 'p0a', 'p1c', 'p1a'), names(theta)))) {
        stop('theta does not contain all of the required parameters:
              lambda0, lambda1, mu0, mu1, q01, q10, p0c, p0a, p1c, p1a')
    }
  } else if(model == 'bd') {
      if(any(!is.element(c('lambda', 'mu'), names(theta)))) {
        stop('theta does not contain all of the required parameters:
             lambda, mu')
      }
  } else if(model == 'classe') { # model depends upon number given parameters
      if(any(!is.element(c('lambda111', 'mu1'), names(theta)))) {
        stop('theta does not contain the minimum required parameters:
             lambda111, mu1. Additional parameters are required for specifying
             more than one trait.')
      }
  } else if(model == 'geosse') {
      if(any(!is.element(c('sA', 'sB', 'sAB', 'xA', 'xB', 'dA', 'dB'),
                         names(theta)))) {
        stop('theta does not contain all of the required parameters:
               sA, sB, sAB, xA, xB, dA, dB')
    }
  } else if(model == 'musse') { # model depends upon number given parameters
      if(any(!is.element(c('lambda1', 'mu1'), names(theta)))) {  
        stop('theta does not contain the minimum required parameters:
             lambda1, mu1. Additional parameters are required for specifying
             more than one trait.')
    }
  } else if(model == 'quasse') {
    if(any(!is.element(c('lambda', 'mu', 'char'), names(theta)))) {
      stop('theta does not contain all of the required parameters:
           lambda, mu, char.')
    }
  } else if(model == 'yule') {
    if(!is.element('lambda', names(theta))) {
      stop('theta does not contain the required parameter: lambda')
    }
  } else {
      stop('model must be set to one of: bisse, bisseness, bd,
           classe, geosse, musse, quasse, yule')
  }
  ## Parse tips information
  tips <- .parse.tips(tips)
  ## Set seed
  if(!is.na(seed)) {
        set.seed(seed)
  }
  ## BiSSE parameter vector
  parms <- unname(theta)
  ## Simulate tree(s)
  result <- trees(parms, type=model, n=nsim, max.taxa=tips$n.tips,
                  max.t=Inf, include.extinct=FALSE, ...)
  return(result)
}
attr(speciation.model, 'name') <- 'speciation.model'
