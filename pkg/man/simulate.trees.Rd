\name{simulate.trees}
\alias{simulate.trees}

\title{
    Simulate Trees from Particle
}
\description{
    \code{simulate.trees} uses the model specified in
    workspace$config to create a list of simulated trees
    that is the length of the number of samples per particle.
    
    Each simulated tree is preprocessed by ladderizing, rescaling,
    and computing the self kernel for the tree.
}

\usage{
    simulate.trees(workspace, theta, seed = NA, ...)
}

\arguments{
    \item{workspace}{
        An object of class 'smc.workspace'.
    }
    \item{theta}{
        A named vector of parameter values.
    }
    \item{seed}{
        Argument to \code{set.seed()}.
    }
    \item{\dots}{
        Additional arguments to pass to config@generator.
    }
}

\details{}

\value{
    A list of objects of class 'phylo', the length of which is
    equal to the number of samples per particle.
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}(5), 1009--1020.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{initialize.smc}} for the use of the output of this function.
}

\examples{
## Set observed tree
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")

## Set SMC configuration and workspace
 config <- list(
        nparticle=10,
        model=const.coalescent,
        nsample=10,
        decay.factor=0.2,
        rbf.variance=2.0,
        sst.control=1.0,
        norm.mode='NONE'
    )

workspace <- init.workspace(t2, config)

## Set Parameters
theta <- c(1)
names(theta) <- c('Ne.tau')

simulate.trees(workspace, theta)
}

\keyword{simulate}