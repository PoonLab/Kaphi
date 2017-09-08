setwd('~/git/Kaphi/')
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
            
                <output spec='NexusOutput' fileName='test.nexus' collapseSingleChildNodes='false'/>
                <output spec='JsonOutput' fileName='test.json'/>
                <output spec='NewickOutput' fileName='test.newick'/>
              </run>
            </beast>"

# <inheritancePostProcessor spec='LineageSampler' reverseTime='false' noClean='false' nSamples='10'>
# </inheritancePostProcessor>
# <inheritancePostProcessor spec="LineageFilter" reactionName="Sampling"/>
# <postSimCondition spec="LeafCountPostSimCondition" nLeaves="{{ ntips }}" exact="false"/>
# <lineageEndCondition spec='LineageEndCondition' nLineages='0'/>


# hash
data <- list(tend = '40',
             beta = '0.001',
             gamma = '0.2',
             phi = '0.01',
             N = (1000-1))

text <- whisker.render(template, data)
write(text, file='pkg/R/temp.xml')



setwd('~/git/MASTER-5.1.1')
system2('java', args=c('-jar MASTER-5.1.1.jar', '../Kaphi/pkg/R/temp.xml'), stdout=F, stderr=F)

require(rjson)
df <- fromJSON(file='test.json')
plot(df$t, df$S, 's', col='green', xlab='Time', ylab='Population size')
points(df$t, df$I, 's', col='red', xlab='Time', ylab='Population size')
points(df$t, df$R, 's', col='blue', xlab='Time', ylab='Population size')

