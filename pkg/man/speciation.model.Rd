\name{speciation.model}
\alias{speciation.model}

\title{
    Library of Standard Speciation Models
}
\description{
    \code{speciation.model} simulates phylogenetic trees for one of eight
    standard speciaiton models (BiSSE, BiSSE-ness, Birth-Death, ClaSSE, GeoSSE,
    MuSSE, QuaSSE, Yule).
    
    Parameter values for the models are given by theta, as set by the user in
    the configuration YAML file.
}

\usage{
    speciation.model(theta, nsim, tips, model, seed = NA)
}

\arguments{
    \item{theta}{
        A vector containing the parameter values for the model.
        See details for the specific parameters required for each model.
    }
    \item{nsim}{
        Number of trees to simulate.
    }
    \item{tips}{
        If integer, the number of tips on the tree. If vector,
        the heights of the tips.
    }
    \item{model}{
        Model under which trees are evolved. Must be one of: 'bisse', 'bisseness',
        'bd', 'classe', 'geosse', 'musse', 'quasse', 'yule'.
    }
    \item{seed}{
        Argument to set seed.
    }
}

\details{
    Binary State Speciation and Extinction (BiSSE) requires 6 parameters: lambda0, 
    lambda1, mu0, mu1, q01, q10. Lambda0 and lambda1 are the rates of specitation in
    character state 0 and 1, respectively. Likewise, mu0 and mu1 are the rates of 
    extinction in each state. q01 is the rate of change from state 0 to 1 and q10 is 
    the rate of change from state 1 to 0.
    
    BiSSE-Node Enhanced State Shift (BiSSE-ness) requires 10 parameters: lambda0, 
    lambda1, mu0, mu1, q01, q10, p0c, p0a, p1c, p1a. In addition to the parameters 
    described for BiSSE, p0c and p1c are the probabilities that a character change 
    from state 0 and 1, respectively, occurs during speciation. p0a and p1a are the 
    probabilities that a state change occuring during specitation is asymmetrical.
    
    Birth-Death (bd) only requires 2 parameters: lambda and mu. All lineages share a 
    fixed probability of speciating or going extinct.
    
    Cladogenetic State change Speciation and Extinction (ClaSSE) requires more parameters
    as the number of states observed increases. For 2 traits the parameters are: lambda111, 
    lambda112, lambda122, lambda211, lambda 212, lambda222, mu1, mu2, q12, q21. The mu 
    and q parameters are comparable to those in the BiSSE model. The first number after 
    lambda indicates the ancestral state of the character. The following 2 number describe 
    the states of the daughter lineages. Thus for lambda111 there is no change in state and
    lambda112 is a state change from 1 to 2 in one of the daughter lineages during 
    speciation. ClaSSE describes cladogenetic change (during speciation) rather than 
    anagenetic change (along a lineage). 
    
    Geographic State Speciation and Extcintion (GeoSSE) requires 7 parameters: sA, sB, sAB,
    xA, xB, dA, dB. sA and sB refer to the rate of speciation in geographic regions A and B.
    sAB decribes the rate of speciation for a lineage that covers both geographic regions.
    xA and xB are the rates of extinction or extirpation from region region A and B. dA and 
    dB are the rate of dispersal/range expansion from region A and B.
    
    Multiple State Speciation and extinction (MuSSE), like ClaSSE, requires parameters that
    depend on the number of states considered. The MuSSE model takes into account the rates 
    of speciation, extinction and state change for the desired number of states. For 3
    states, the parameters required are: lambda1, lambda2, lambda3, mu1, mu2, mu3, q12, q13,
    q21, q23, q31, q32.
    
    Quantitative state Speciation and Extinction requires 3 parameters: lambda, mu, and char.
    char is a model of character evolution. For this model, lambda and mu must be functions
    of the character state, x. Thus the rates of speciation and extinction are lambda(x) and
    mu(x).
    
    Yule model only requires one parameter: lambda, the rate of speciation. 
}

\value{
    An object of class phylo or a list of phylo trees.
}

\references{
BiSSE
	Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary character's effect on speciation and extinction. Systematic Biology 56: 701-710.

	FitzJohn, R. G., Maddison, W. P., & Otto, S. P. (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic biology, 58(6), 595-611.

BiSSE-ness
	Magnuson-Ford K. and Otto S.P. 2012. Linking the investigations of character evolution and species di-versification. The American Naturalist

Birth-Death
	Harmon, L. J., Weir, J. T., Brock, C. D., Glor, R. E., & Challenger, W. (2008). GEIGER: investigating evolutionary radiations. Bioinformatics, 24(1), 129-131.

ClaSSE
	Goldberg, E. E., & Igić, B. (2012). Tempo and mode in plant breeding system evolution. Evolution, 66(12), 3701-3709.

GeoSSE
	Goldberg, E. E., Lancaster, L. T., & Ree, R. H. (2011). Phylogenetic inference of reciprocal effects between geographic range evolution and diversification. Systematic Biology, 60(4), 451-465.

MuSSE
	FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6), 1084-1092.

QuaSSE
	FitzJohn, R. G. (2010). Quantitative traits and diversification. Systematic biology, syq053.

Yule
  Steel, M., & McKenzie, A. (2001). Properties of phylogenetic trees generated by Yule-type speciation models. Mathematical biosciences, 170(1), 91-112.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{trees}} for tree simulating function from Diversitree package.
}

\examples{
## Example use of function for BiSSE model tree simulation. Each model requires
## differenet parameters to be set in the configuration YAML file (see details).

require(diversitree) # If not loaded already
set.seed(1)

## If running from R GUI:
setwd('~/git/Kaphi')

## Load configuration file:
config <- load.config('pkg/examples/example-bisse.yaml')
config <- set.model(config, 'bisse')

## Sample Priors for parameter values
theta <- sample.priors(config)

## Simulate trees:
sim.trees <- speciation.model(theta, nsim=5, tips=20, 'bisse', seed=NA)
}

\keyword{bisse}
\keyword{bisseness}
\keyword{bd}
\keyword{classe}
\keyword{geosse}
\keyword{musse}
\keyword{quasse}
\keyword{yule}
\keyword{speciation}
\keyword{model}
