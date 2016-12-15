# Kaphi

## Brief description
* [Phylodynamic inference](https://en.wikipedia.org/wiki/Viral_phylodynamics) is the fitting of models to the shape of a [phylogenetic tree](https://en.wikipedia.org/wiki/Phylogenetic_tree) in order to reconstruct the historical processes that produced the tree.
* Kaphi is an [R package](https://cran.r-project.org/) for fitting models to the shapes of phylogenetic trees.  
* It uses [approximate Bayesian computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation) as a likelihood-free approach to fit models.
* The shapes of simulated phylogenies are compared to the data with a [kernel method](https://en.wikipedia.org/wiki/Kernel_method).


## Phylodynamic inference

A [phylogeny](https://en.wikipedia.org/wiki/Phylogenetic_tree) is a tree-based hypothesis about how different populations are related by common ancestors.
The idea behind phylodynamics is that the 


## Approximate Bayesian computation
Approximate Bayesian computation (ABC) is a class of methods in statistical computing.  
The basic idea behind ABC is that we can fit a model by simulating data sets under different model parameter settings, comparing these simulations to the actual data, and then keeping the parameter values that yield the most realistic simulations. 
The key advantage of ABC over conventional likelihood-based methods for model fitting is that we don't have to calculate the likelihood of the model.  
Depending on the complexity of the model, calculating its likelihood may be time-consuming or simply not feasible.
On the other hand, it is usually much easier to use the model to simulate data sets. 

ABC requires a method for sampling parameter values from its approximation of the posterior distribution.
Sequential Monte Carlo (SMC) or particle filtering is a method in which a population of particles (that represent vectors of model parameter values) is initialized from a prior distribution, and then iteratively updated so they converge to the posterior distribution. 



## Okay, what *is* Kaphi?

Kaphi is an R package 


## What's with the name?
Kaphi is an abbreviation of "**K**ernel-embedded **A**BC-SMC for **ph**ylodynamic **i**nference". 
An earlier version that used ABC-MCMC was called *Kamphir*, but we've dropped MCMC in favour of SMC. 
Besides, my partner didn't like the old name because it reminded her of [mothballs](https://en.wikipedia.org/wiki/Camphor).  
Fun fact: apparently it's also Sanskrit for "coffee". 


## Financial support
This work is funded by the Government of Canada through [Genome Canada](https://www.genomecanada.ca/) and the [Ontario Genomics Institute](http://www.ontariogenomics.ca/) (OGI-131).