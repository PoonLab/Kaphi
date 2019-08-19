require(Kaphi)
require(reticulate)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# bring in python for python scripts
use_python("~/miniconda2/bin/python2.7", required = T)
# py_config()

# python scripts
source_python('~/Documents/coevolutionMisc/nestedCoalescent.py')
source_python('~/Documents/coevolutionMisc/phyloK2.py')

# host tree file name
host.tree <- '/home/tng92/Documents/coevolutionMisc/true-trees/Hostsimulatedtree.nwk'

# load configuration file
config <- load.config('~/git/Kaphi/pkg/examples/example-coevolution.yaml')
config <- set.model(config, 'coevolution')

# simulate target tree
theta <- c(L=0.1, M=0, P=0.01)  # this is the true value
starting.trees <- main(host.tree, 0.1, 0, 0.01, FALSE)
host.tree <- read.tree(text=starting.trees[[1]])
obs.tree <- read.tree(text=starting.trees[[2]])
obs.tree <- parse.input.tree(obs.tree, config)

## now let's estimate that posterior distribution!

# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
res <- run.smc(ws, trace.file='pkg/examples/example-coevolution.tsv', nthreads = 10, model='coevolution', verbose=TRUE)

