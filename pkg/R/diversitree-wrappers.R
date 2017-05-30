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
speciation.model <- function(theta, nsim, tips, model, seed=NA) {
    "
    theta : a vector containing the parameter values for the model
        @param lambda(x) : Rate of speciation (in state x)
        @param mu(x) : Rate of extinction (in state x)
        @param qxy : Rate of change from state x to y
        @param pxc : Probabilitiy that a character change from ancetral
                     state x occures during speciation
        @param pxa : Probability that a character change from ancestral
                     state x during speciation is asymetrical
        @param
    nsim : number of trees to simulate
    tips : if integer, the number of tips; if vector, the heights of the tips
    model : model under which trees are evolved. One of:
              - bisse       -- Binary State Speciation and Extinction
              - bisseness   -- BiSSE-Node Enhanced State Shift
              - bd          -- Birth-Death
              - classe      -- Cladogenetic State change Speciation and Extinction
              - geosse      -- Geographic State Speciation and Extinction
              - musse       -- Multi-state Speciation and Extinction
              - quasse      -- Quantitative state Specatiation and Extinction
              - yule        -- Yule model

    Use diversitree to simulate trees under one of the 'SSE' models or a simple
    character indepdenent birth-death model.
    "
  ## Validate arguments
  #if (model == bisse) {
  #  if(all(!is.element(c('lambda0', 'lambda1', 'mu0', 'mu1', 'q01', 'q10'), names(theta)))) {
  #    stop('theta does not contain required parameters')
  #  }
  #}  
  if(!is.element(model, c('bisse', 'bisseness', 'bd', 'classe', 'geosse', 'musse', 'quasse', 'yule'))) {
    stop('model must be set to one of: bisse, bisseness, bd, classe, geosse, musse, quasse, yule')
  }
  ## Set seed
  if(!is.na(seed)) {
        set.seed(seed)
  }
  ## BiSSE parameter vector
  parms <- unname(theta)
  ## Simulate tree(s)
  result <- trees(parms, type=model, n=nsim, max.taxa=tips, 
                  max.t=Inf, include.extinct=FALSE)
  return(result)  
}

set.model.name <- function(speciation.model) {
  "
  Use of this function depends on whether the model is going to be
  called 'speciation.model' or by the specific model set by the user.
  
  If called 'speciation.model', `attr(speciation.model, 'name') <- 'speciation'`
  will suffice.
  If called by the specific model used, this function can serve to set the name
  attribute of the model to the correct model name
  "
}
## Figure out how to handle these:
attr(bisse, 'name') <- "bisse"
attr(bisseness, 'name') <- "bisseness"
## ...

