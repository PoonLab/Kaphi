\name{initialize.smc}
\alias{initialize.smc}

\title{
    Initialize Sequential Monte Carlo Workspace
}

\description{
    \code{initialize.smc} updates a 'smc.workspace' object with
    the data necessary for use by \code{run.smc}.
    
    This involves populating the contents of variables
    workspace$particles, workspace$weights, workspace$sim.trees,
    and workspace$dists with the initial values for \code{run.smc}.
}

\usage{
    initialize.smc(ws, ...)
}


\arguments{
    \item{ws}{
        An object of class 'smc.workspace'.
    }
    \item{\dots}{
        Additional arguments to pass to config@generator via simulate.trees().
    }
}

\details{

}

\value{
    An object of class 'smc.workspace' with updated contents for
    the variables: particles, weights, sim.trees, and dists.
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}(5), 1009--1020.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{simulate.trees}} for the creation of sim.trees.
    \code{\link{distance}} for the computation of distances in dists.
}

\examples{
## Set observed tree
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")

## Set SMC configuration (will need to setwd() if running in GUI)
config <- load.config('tests/fixtures/test-coalescent.yaml')

## Set model
config <- set.model(config, const.coalescent)

## Initialize workspace
ws <- init.workspace(t2, config)

ws <- initialize.smc(ws)
}

\keyword{initialize}
\keyword{ ~kwd2 }