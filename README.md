# Kaphi

## Brief description
* [Phylodynamic inference](https://en.wikipedia.org/wiki/Viral_phylodynamics) is the fitting of models to the shape of a [phylogenetic tree](https://en.wikipedia.org/wiki/Phylogenetic_tree) in order to reconstruct the historical processes that produced the tree.
* Kaphi is an [R package](https://cran.r-project.org/) for fitting models to the shapes of phylogenetic trees.  
* It uses [approximate Bayesian computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation) as a likelihood-free approach to fit models.
* The shapes of simulated phylogenies are compared to the data with a [kernel method](https://en.wikipedia.org/wiki/Kernel_method).


## Approximate Bayesian computation
Approximate Bayesian computation (ABC) is a class of methods in statistical computing that generally fit a model to the data by simulating data sets under different versions of the model, comparing these simulations to the actual data, and then keeping the parameter values that yield the most realistic simulations.

## What's with the name?
Kaphi is an abbreviation of "**K**ernel-embedded **A**BC-SMC for **ph**ylodynamic **i**nference".  An earlier version that used ABC-MCMC was called *Kamphir*, but we've dropped MCMC in favour of SMC - besides, my partner didn't like the old name because it reminded her of [mothballs](https://en.wikipedia.org/wiki/Camphor). 

## Financial support
This work is funded by the Government of Canada through [Genome Canada](https://www.genomecanada.ca/) and the [Ontario Genomics Institute](http://www.ontariogenomics.ca/) (OGI-131).
