# Kaphi

## Brief description
* [Phylodynamic inference](https://en.wikipedia.org/wiki/Viral_phylodynamics) is the fitting of models to the shape of a [phylogenetic tree](https://en.wikipedia.org/wiki/Phylogenetic_tree) in order to reconstruct the historical processes that produced the tree.
* Kaphi is an [R package](https://cran.r-project.org/) for fitting models to the shapes of phylogenetic trees.  
* It uses [approximate Bayesian computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation) as a likelihood-free approach to fit models.
* The shapes of simulated phylogenies are compared to the data with a [kernel method](https://en.wikipedia.org/wiki/Kernel_method).


## Installation

Instructions on installing Kaphi are provided in [INSTALL.md](INSTALL.md).


## License

Kaphi is released under the [GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.en.html).


## Phylodynamic inference

A [phylogeny](https://en.wikipedia.org/wiki/Phylogenetic_tree) is a tree-based hypothesis about how different populations are related by common ancestors.  The general idea behind phylodynamics is that this tree has been shaped by a history of biological processes, such as speciation or epidemiology (when the tree relates cases of an infectious virus or bacteria).  The objective of phylodynamics is to reconstruct this biological process by fitting a model of that process to the shape of a tree.  


## Approximate Bayesian computation
Approximate Bayesian computation (ABC) is a class of methods in statistical computing.  The basic idea behind ABC is that we can fit a model by simulating data sets under different model parameter settings, comparing these simulations to the actual data, and then keeping the parameter values that yield the most realistic simulations.  The key advantage of ABC over conventional likelihood-based methods for model fitting is that we don't have to calculate the likelihood of the model.  Depending on the complexity of the model, calculating its likelihood may be time-consuming or simply not feasible.  On the other hand, it is usually much easier to use the model to simulate data sets. 

The trick to using ABC is to have a measure of how similar a simulation is to the actual data.  This is not a trivial problem because the shape of a tree is a complex property that is difficult to reduce down to a number or set of numbers in a meaningful way.  For this purpose, Kaphi uses a technique from machine learning known as the [kernel method](https://en.wikipedia.org/wiki/Kernel_method).  The particular kernel method used in Kaphi breaks trees down into fragments and then basically counts the number of similar-shaped fragments that two trees have in common.


## Okay, what does Kaphi do?

Kaphi is an R package for fitting models to tree shapes using a combination of kernel methods and ABC.  It is essentially a re-implementation and expansion of a set of C programs written by @rmcclosk called [netabc](https://github.com/rmcclosk/netabc).  This procedure can be used to fit conceivably *any* model to the shape of a tree, as long as that model can be used to simulate trees.  Kaphi is designed to fit any model that can be expressed as systems of ordinary differential equations, *i.e.*, compartmental 'mass action' models.  Since Kaphi is derived from *netabc*, it can also fit models that generate trees as an epidemic process on a graph (network). 


## What's with the name?
Kaphi is an abbreviation of "**K**ernel-embedded **A**BC-SMC for **ph**ylodynamic **i**nference".  An earlier version that used ABC-MCMC was called *Kamphir*, but we've dropped MCMC in favour of SMC.  Besides, my partner didn't like the old name because it reminded her of [mothballs](https://en.wikipedia.org/wiki/Camphor).  Fun fact: apparently "Kaphi" is also Sanskrit for "coffee". 


## Financial support
This work is funded by the Government of Canada through [Genome Canada](https://www.genomecanada.ca/) and the [Ontario Genomics Institute](http://www.ontariogenomics.ca/) (OGI-131).
