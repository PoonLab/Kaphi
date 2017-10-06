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
  
  res.trees <- .call.master(theta=theta, nsim=nsim, tips=tips, seed=seed, tsample=tsample)
      
  return(res.trees)
}
attr(epidem.model, 'name') <- "epidem.model"  # satisfies requirement in smcConfig.R set.model() function



.call.master <- function(theta, nsim, tips, seed=NA, tsample) {
  require(whisker, quietly=TRUE)
  tree.file <- tempfile(pattern='temp', fileext='.newick') 
  
  ## hash
  data <- list(tend = as.character(theta$t.end),
               beta = as.character(theta$beta),
               gamma = as.character(theta$gamma),
               phi = as.character(theta$phi),
               N = as.character(theta$N - 1),
               #seed = as.character(NA),
               nsample = as.character(nsim),
               ntips = as.character(tips),
               tsampl = as.character(tsample),
               tempNewick = tree.file         
               )
  # current problem is tempfile() is relative path (ie. '/tmp/RtmpNhrcQi/temp34e847bc77fa.newick')
  # XML template will not render with the '/', needs only the fileName (ie. # current problem is tempfile() is relative path (ie. 'temp34e847bc77fa.newick')
  
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

                <lineageEndCondition spec='LineageEndCondition' population='@I' nLineages='0' isRejection='true'>
                </lineageEndCondition>
               
                <inheritancePostProcessor spec='LineageSampler' samplingTime='{{ tsampl }}'>
                  <populationSize spec='PopulationSize' population='@I_sample' size='{{ ntips }}'/>
                </inheritancePostProcessor>
            
                <output spec='NewickOutput' fileName='{{ tempNewick }}' collapseSingleChildNodes='true'/>
              </run>
            </beast>"
  
  
  #  <populationEndCondition spec='PopulationEndCondition' threshold='{{ ntips }}' exceedCondition='true' population='@I_sample'/>
  
  
  
  ## generate temporary XML
  text <- whisker.render(template, data)
  tempname <- tempfile(pattern='temp', fileext= '.xml')
  write(text, file=tempname)
   
  ## system call to MASTER with temporarily generated XML
  system2('java', args=c('-jar ../MASTER-5.1.1/MASTER-5.1.1.jar', paste0(tempname)), stdout=F) #, stdout=F, stderr=F
  

  if (!file.exists(tree.file)) {
    # should re-simulate trees
    return(0)
  }
  ## read Newick, reset to Kaphi directory, and send tree back to user
  # casting result as a multiPhylo object
  trees <- tryCatch({
    read.tree(file=tree.file, keep.multi=TRUE)
    }, warning = function(w) {
      cat('The warning ', w)
      # should re-simulate trees, do not continue feeding the same theta params over and over
      return(0)
    }, error = function(e) {
      cat('The error: ', print(paste0(e)))
      # should re-simulate trees, do not continue feeding the same theta params over and over
      return(0)
    }, finally = {
      unlink(c(tempname, tree.file))
    })
  
  if (length(trees) != nsim) {
    # if "Warning: Newick writer skipping empty graph" --> one of newick trees eliminated
    # throws "Error in x[[i]]: subscript out of bounds."
    cat('Less than 5 trees retained. Re-simulating..')
    return(0)
  }
  
  return(trees)
}

