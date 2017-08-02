***(work in progress)***

# Configuration File Format

Kaphi settings are read in from a user-created [YAML](http://yaml.org) file. This document outlines the format that is required by the YAML parser.

The examples given within this document are specific to the birth-death speciation model.
Further examples of YAML files can be found in *Kaphi/pkg/examples*.

### Prior Distributions
A prior distribution is chosen for each model parameter, from which the particle values are sampled. 

Each parameter's prior is specified by the following:
* *parameter name*
* *dist* - one of the [distributions in the stats package](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html), without a prefix. 
* *hyperparameters* - a list of the hyperparameter names and associated values for the chosen distribution.

The example below is assigning a gamma distribution prior to the parameter lambda and a uniform distribution to mu.

```YAML
priors:
  lambda:
    dist: 'gamma'  
    hyperparameters:
    - rate: 1.0  
    - shape: 2.0
  mu:
    dist: 'unif'
    hyperparameters:
    - min: 0.0
    - max: 0.05
constraints:
  - 'mu < lambda'
```
#### Constraints
The priors are followed by the constraints section. 
This section allows the user to enter a string containing any conditions that must be upheld during sampling (or will otherwise leave vacant).

As shown in the example above, constraints are of the form: `"param1 operator param2"`. 

### Proposal Distributions
A proposal distribution is specified for each parameter. The values sampled from these distributions dictate how a particle is updated.

The format is similar to that of the priors, with the exception that *parameters* are specified rather than *hyperparameters*.
```YAML
proposals:
  lambda:
    dist: 'norm'
    parameters:
    - mean: 0
    - sd: 0.1
  mu:
    dist: 'norm'
    parameters:
    - mean: 0
    - sd: 0.1
```

### SMC Settings

The Sequential Monte Carlo (SMC) settings allow for adjustments to be made to different functions within smcABC.R.
```YAML
smc:
  nparticle: 10
  nsample: 5
  ess.tolerance: 50.0
  final.accept.rate: 0.01
  final.epsilon: 0.01
  quality: 0.95
  step.tolerance: 1.e-4
```

### Distance Metric Settings
Distance metrics are used to compare two trees. 
Kaphi was created to be used with the kernel distance, but many other distaces 

```YAML
distances:
  'kernel.dist':
    package: 'Kaphi'
    weight: 0.8
    decay.factor: 0.2
    rbf.variance: 100.0
    sst.control: 1.0
  'sackin':
    package: 'Kaphi'
    weight: 0.05      # making up the rest of these parameters, just for testing purposes
    use.branch.lengths: FALSE
  'colless':
    package: 'Kaphi'
    weight: 0.10
  'gammaStat':
    package: 'ape'
    weight: 0.05
```

```YAML
distances:
  '0.8*kernel.dist(x,y,decay.factor=0.2,rbf.variance=100.0,sst.control=1.0)+0.1*sackin+0.3*colless':
```


#### Valid Distance Metrics
| Tree Statistic | Package | Description |
|----------------|---------|-------------|
| tree.kernel*| Kaphi |The kernel distance between two trees.|
| ~~nLTT~~| ~~Kaphi~~|*Does not currently work*|
| sackin | Kaphi |Sackin index.|
| colless | Kaphi |Colless imbalance number.|
| ~~cophenetic~~| ~~Kaphi~~ |*Does not currently work*|
| ladder.length | Kaphi |Max ladder length.|
| IL.nodes | Kaphi |Number of internal nodes with one leaf.|
| tree.width | Kaphi |Max width divided by max depth.|
| max.delta.width | Kaphi |Max difference in in width between two levels.|
| n.cherries | Kaphi |Number of cherries(node with two leaves).|
| prop.unbalanced | Kaphi |Proportion of unbalanced subtrees.|
| avg.unbalance | Kaphi |Average ratio of unbalanced subtrees.|
| pybus.gamma | Kaphi |Pybus' gamma statistic|
| internal.terminal.ratio | Kaphi |Ratio of internal to terminal branches.|
| balance** | ape |For each node of the tree, the numbers of tips on each of its daughter-branches.|
| cophenetic.phylo** | ape |Pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.|
| dist.nodes** | ape |Like `cophenetic.phylo`, but includes internal nodes.|
| dist.topo* | ape |Topological distance between two trees using the method from Penny & Hendy (1985).|
| ~~gammaStat~~|~~ape~~|(pybus.gamma in Kaphi)|
| avgladder | phyloTop |Mean size of ladders in the tree.|
| ~~cherries~~ | ~~yloTop~~|(n.cherries in Kaphi)|
| ~~colless.phylo~~ | ~~phyloTop~~|(colless in Kaphi)|
| getDepths** | phyloTop |Returns a list of two vectors: `tipDepths` and `nodeDepths`.|
| ~~ILnumber~~ | ~~phyloTop~~ |(IL.nodes in Kaphi)|
| ~~maxHeight~~ | ~~phyloTop~~ |(ladder.length in Kaphi)|
| pitchforks | phyloTop |Number of clades with three tips. |
| ~~sackin.phylo~~ | ~~phyloTop~~ |(sackin in Kaphi)|

\* requires two trees.
** output is non-scalar.

## Template
```YAML
priors:
  <parameter_1>:
    dist: ''
    hyperparameters:
    -

constraints: 
 - ''
 
proposals:
  <parameter_1>:
    dist: ''
    parameters:
    - 

smc:
  nparticle: 100
  nsample: 5
  ess.tolerance: 50.0
  final.accept.rate: 0.05
  final.epsilon: 0.05
  quality: 0.95
  step.tolerance: 1.e-5
  norm.mode: 'NONE'

distances:
  'kernel.dist':
    package: 'Kaphi'
    weight: 0.8
    decay.factor: 0.2
    rbf.variance: 100.0
    sst.control: 1.0
```