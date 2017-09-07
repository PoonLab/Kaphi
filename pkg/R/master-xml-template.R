require(whisker)

# template
template <- 
           "<beast version='2.0' namespace= 'master
                                             :master.model
                                             :master.conditions
                                             :master.outputs
                                             :master.postprocessors'>
              <run spec='InheritanceTrajectory'
                   samplePopulationSizes='true'
                   verbosity='2'
                   simulationTime='10'
                   seed='42'>
            
                <model spec='Model' id='model'>
                  <population spec='Population' id='V' populationName='V' />
            
                  <reaction spec='Reaction' reactionName='birth' rate='1.0'>
                    {{{reaction}}}
                  </reaction>
                </model>
            
                <initialState spec='InitState'>
                  <lineageSeed spec='Individual' population='@V' />
                </initialState>
            
                <inheritancePostProcessor spec='LineageSampler' reverseTime='false' noClean='false' nSamples='10'>
            <!--
                  <populationSize spec='PopulationSize' size='10'>
                    <population spec='Population' populationName='T'/>
                  </populationSize>
             -->
            <!--
                  <populationSize spec='PopulationSize' size='10'>
                    <population spec='Population' populationName='V'/>
                  </populationSize>
             -->
                </inheritancePostProcessor>
            
                <output spec='NexusOutput' fileName='test.nexus' collapseSingleChildNodes='false' />
                <output spec='JsonOutput' fileName='test.json' />
              </run>
            </beast>"

# hash
data <- list(reaction = 'V -> V + V')


text <- whisker.render(template, data)

setwd('~/git/Kaphi/')
write(text, file='pkg/R/temp.xml')
