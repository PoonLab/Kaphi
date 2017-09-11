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


epidem.model <- function(theta, nsim, tips, model='sir.nondynamic', seed=NA, labels=NA) {
  theta <- as.list(theta)
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  trees <- replicate(nsim, # num of simulations
                     .call.master(theta),
                     simplify=FALSE
                     )
  # cast result as a multiPhylo object
  class(trees) <- "multiPhylo"
  return(trees)
}
attr(compartmental.model, 'name') <- "epidem.model"  # satisfies requirement in smcConfig.R set.model() function



.call.master <- function(theta) {
  setwd('~/git/Kaphi/pkg/R')
  require(whisker)

  ## hash
  data <- list(tend = as.character(theta$t.end),
               beta = as.character(theta$beta),
               gamma = as.character(theta$gamma),
               phi = as.character(0.0002),
               N = as.character(theta$N - 1))
  
  ## XML template
  template <- 
             "<beast version='2.0' namespace= 'master
                                               :master.model
                                               :master.conditions
                                               :master.outputs
                                               :master.postprocessors'>
                <run spec='InheritanceTrajectory'
                     samplePopulationSizes='true'
                     verbosity='2'
                     simulationTime='{{ tend }}'
                     seed='42'>
              
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
                      I -> I_sample
                    </reaction>
                  </model>
              
                  <initialState spec='InitState'>
                    <populationSize spec='PopulationSize' population='@S' size='{{ N }}'/>
                    <lineageSeed spec='Individual' population='@I'/>
                  </initialState>
              
                  <output spec='JsonOutput' fileName='temp.json'/>
                  <output spec='NewickOutput' fileName='temp.newick'/>
                </run>
              </beast>"
  
  # <inheritancePostProcessor spec='LineageSampler' reverseTime='false' noClean='false' nSamples='10'>
  # </inheritancePostProcessor>
  # <inheritancePostProcessor spec="LineageFilter" reactionName="Sampling"/>
  # <postSimCondition spec="LeafCountPostSimCondition" nLeaves="{{ ntips }}" exact="false"/>
  # <lineageEndCondition spec='LineageEndCondition' nLineages='0'/>
  
  ## generate temporary XML
  text <- whisker.render(template, data)
  write(text, file='temp.xml')
  
  ## system call to MASTER with temporarily generated XML
  system2('java', args=c('-jar ../../../MASTER-5.1.1/MASTER-5.1.1.jar', 'temp.xml'), stdout=F, stderr=F)
  
  ## read Newick, reset to Kaphi directory, and send tree back to user
  tree <- read.tree(file='temp.newick')
  setwd('~/git/Kaphi/')
  tree
}



## plotting simulation from json file
require(rjson)
df <- fromJSON(file='test.json')
plot(df$t, df$S, 's', col='green', xlab='Time', ylab='Population size')
points(df$t, df$I, 's', col='red', xlab='Time', ylab='Population size')
points(df$t, df$R, 's', col='blue', xlab='Time', ylab='Population size')

