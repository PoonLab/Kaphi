require(Kaphi)

# I need to set this if running from R GUI
setwd('~/git/Kaphi')

# load configuration file (assumes R was launched from Kaphi root dir)
#config <- load.config('pkg/examples/example-bisse.yaml')
config <- load.config('pkg/examples/phyloK3_2-example.yaml')
config <- set.model(config, 'bisse')

# simulate target tree
theta <- c(lambda0=0.2, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)  # this is the true value
set.seed(50)
obs.tree <- speciation.model(theta, nsim=1, tips=50, model='bisse')[[1]]
obs.tree <- parse.input.tree(obs.tree, config)
write(write.tree(obs.tree),"~/git/Coevolution/examples/obstree_lambda0.nwk")

# Visualize parameter identifiability

# calculate kernel distances for varying lambda0
png(filename = "~/Documents/BiSSE/kernel50_1_Lambda0.png",width=900,height=900,res=120)
x <- runif(1000, 0, 1)
v <- x
res <- sapply(x, function(val) {
  theta <- c(lambda0=val, lambda1=0.1, mu0=0.003, mu1=0.003, q01=0.03, q10=0.01)
  sim.trees <- speciation.model(theta, nsim=1, tips=50, model='bisse')
})

for (i in res){
  write(write.tree(i), "~/git/Coevolution/examples/python_lambda0.nwk", append = TRUE)}