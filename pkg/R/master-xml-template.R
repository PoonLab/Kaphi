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


epidem.model <- function(theta, nsim, tips, tsample, model='epidemic', seed=NA, labels=NA) {
  th.args <- names(theta)
  if (length(th.args) < 5 || any(!is.element(c('t.end', 'N', 'beta', 'gamma', 'phi'), th.args))) {
    stop("'theta' does not hold Kaphi-compatible parameters")
  }
  theta <- as.list(theta)
  
  # check if there are tip labels (dates and states) to incorporate
  ntips <- .parse.tips(tips)                   # function is written in smcConfig.R
  if (!is.na(seed)) { set.seed(seed) }
  
  trees <- .call.master(theta, nsim=nsim, tips=tips, seed=seed, tsample=tsample)
      
  return(trees)
}
attr(epidem.model, 'name') <- "epidem.model"  # satisfies requirement in smcConfig.R set.model() function



.call.master <- function(theta, nsim, tips, seed=NA, tsample) {
  setwd('~/git/Kaphi/pkg/R')
  require(whisker, quietly=TRUE)

  ## hash
  data <- list(tend = as.character(theta$t.end),
               beta = as.character(theta$beta),
               gamma = as.character(theta$gamma),
               phi = as.character(theta$phi),
               N = as.character(theta$N - 1),
               #seed = as.character(NA),
               nsample = as.character(nsim),
               ntips = as.character(tips),
               tsampl = as.character(tsample))
  
  ## XML template
  template <- 
            "<beast version='2.0' namespace= 'master
                                             :master.model
                                             :master.conditions
                                             :master.outputs
                                             :master.postprocessors'>
              <run spec='InheritanceEnsemble' 
                   simulationTime='{{ tend }}'
                   samplePopulationSizes='true' 
                   nTraj='{{ nsample }}'
                   verbosity='1'>
            
                <model spec='Model' id='model'>
                  <population spec='Population' id='S' populationName='S' />
                  <population spec='Population' id='I' populationName='I' />
                  <population spec='Population' id='R' populationName='R' />
                  <population spec='Population' id='I_sample' populationName='I_sample' />
            
                  <reaction spec='Reaction' reactionName='Infection' rate='{{ beta }}'>
                    S + I -> 2I
                  </reaction>
                  <reaction spec='Reaction' reactionName='Recovery' rate='{{ gamma }}'>
                    I -> R
                  </reaction>
                  <reaction spec='Reaction' reactionName='Sampling' rate='{{ phi }}'>
                    I:1 -> I_sample:1
                  </reaction>
                </model>
            
                <initialState spec='InitState'>
                  <populationSize spec='PopulationSize' population='@S' size='{{ N }}'/>
                  <populationSize spec='PopulationSize' population='@I_sample' size='0'/>
                  <populationSize spec='PopulationSize' population='@R' size='0'/>
                  <lineageSeed spec='Individual' population='@I'/>
                </initialState>

                <lineageEndCondition spec='LineageEndCondition' population='@I' nLineages='0' isRejection='true'/>
               
                <inheritancePostProcessor spec='LineageSampler' samplingTime='{{ tsampl }}'>
                  <populationSize spec='PopulationSize' population='@I_sample' size='{{ ntips }}'/>
                </inheritancePostProcessor>
            
                <output spec='NewickOutput' fileName='temp.newick' collapseSingleChildNodes='true'/>
              </run>
            </beast>"
  
  
  #  <populationEndCondition spec='PopulationEndCondition' threshold='{{ ntips }}' exceedCondition='true' population='@I_sample'/>
  
  
  
  ## generate temporary XML
  text <- whisker.render(template, data)
  write(text, file='temp.xml')
  
  ## system call to MASTER with temporarily generated XML
  system2('java', args=c('-jar ../../../MASTER-5.1.1/MASTER-5.1.1.jar', 'temp.xml'), stdout=F, stderr=F)
  
  ## read Newick, reset to Kaphi directory, and send tree back to user
  # casting result as a multiPhylo object
  trees <- read.tree(file='temp.newick', keep.multi=TRUE)
  setwd('~/git/Kaphi/')
  trees
}

