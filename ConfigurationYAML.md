***(work in progress)***

# Configuration File Format

Kaphi settings are read in from a user-created [YAML](http://yaml.org) file. This document outlines the format that is required by the YAML parser.

The examples given within this document are specific to the birth-death speciation model.
Additional examples of YAML files can be found in *Kaphi/pkg/examples*.

### General Template
```YAML
priors:
  parameter_1:
    dist: ''
    hyperparameters:
    -

constraints: 
 - ''
 
proposals:
  parameter_1:
    dist: ''
    parameters:
    - 

smc:
  nparticle: 1000
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

### Prior Distributions
A [prior distribution](https://en.wikipedia.org/wiki/Prior_probability) is chosen for each model parameter, from which the particle values are sampled. 

Each parameter's prior is specified by the following:
* *parameter name*
* *dist* - one of the [distributions in the stats package](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html). Only the root of the distribution should be used (*ex.* for the gamma distribution, write 'gamma' instead of 'pgamma', 'rgamma', 'dgamma', 'qgamma').
* *hyperparameters* - a list of the hyperparameter names and associated values for the chosen distribution.

The example below is assigning a gamma distribution prior to the parameter lambda and a uniform distribution to mu.

#### Distributions in the Stats package
| Distribution       | Keyword  | Parameters          |
|--------------------|----------|---------------------|
| Beta               | beta     | shape1, shape2, ncp |
| Binomial           | binom    | size, prob          |
| Cauchy             | cauchy   | location, scale     |
| Chi-squared        | chisq    | df, ncp             |
| Exponential        | exp      | rate                |
| F distribution     | f        | df1, df2, ncp       |
| Gamma              | gamma    | rate, shape         |
| Geometric          | geom     | prob                |
| Hypergeometic      | hyper    | m, n, k             |
| Log-normal         | lnorm    | meanlog, sdlog      |
| Multinomial        | multinom | size, prob          |
| Negative bionomial | nbinom   | size, prob, mu      |
| Normal             | norm     | mean, sd            | 
| Poisson            | pois     | lambda              |
| Student's t        | t        | df, ncp             |
| Uniform            | unif     | min, max            |
| Weibull            | weibull  | shape, scale        |

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
A [proposal distribution](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm) is specified for each parameter. The values sampled from these distributions dictate how a particle is updated.

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
  nparticle: 1000
  nsample: 5
  ess.tolerance: 50.0
  final.accept.rate: 0.05
  final.epsilon: 0.05
  quality: 0.95
  step.tolerance: 1.e-5
  norm.mode: 'NONE'
```
| Parameter         | Description                                                           |
|-------------------|-----------------------------------------------------------------------|
| nparticle         | The number of particles to be run.                                    |
| nsample           | The number of samples per particle.                                   |
| ess.tolerance     | [Effective Sample Size](https://www.johndcook.com/blog/2017/06/27/effective-sample-size-for-mcmc/) tolerance. |
| final.accept.rate | The accept rate that will stop the algorithm.                         |
| final.epsilon     | The value of epsilon that will stop the algorithm.                    |
| quality           | A tuning parameter (alpha) that corresponds to the ratio of current effective sample size (ESS) over previous ESS.  If alpha is near 1, then ESS is held constant over time and particles stay diverse but are slower to converge.  If alpha is near 0, then particles can collapse to a single point that is not necessarily accurate. |
| step.tolerance    | The convergence tolerance for finding the next epsilon.               |

### Distance Metric Settings
Distance metrics are used to compare two trees. 
Kaphi was created to be used with the kernel distance, but many other distances are accepted (see [Distance Metrics](https://github.com/PoonLab/Kaphi/blob/master/DistanceMetrics.md)). Many difference distance metrics may be used in conjunction, as shown in the examples below.

There are two ways to specify distance metrics: as a YAML dictionary (first example) or as a string indicating a linear combination of distances. (second example).

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
    weight: 0.05
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
  '0.8*kernel.dist(x,y,decay.factor=0.2,rbf.variance=100.0,sst.control=1.0)+0.1*sackin+0.3*colless'
```
