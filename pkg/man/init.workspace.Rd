\name{init.workspace}
\alias{init.workspace}
%\alias{parse.input.tree}
%\alias{parse.labels}

\title{
    Initialize the SMC Workspace
}

\description{
    \code{init.workspace} initializes a workspace to be used
    by \code{run.smc} and related functions in smcABC.R.
    
    The function packages the observed tree, SMC configurations
    (from smcConfig.R) and derived data into a 'smc.workspace'
    object.
}

\usage{
    init.workspace(obs.tree, config, regex = NA)
}

\arguments{
    \item{obs.tree}{
        An object of class "phylo".
    }
    \item{config}{
        An object of class "smc.config".
    }
    \item{regex}{
        A character string containing a \link{regular expression}
        to be used by \code{parse.labels}. If NA, tip labels are
        listed as "NA".
    }
}

\details{}

\value{
    A list of class "smc.workspace" containing:
    \item{obs.tree}{
        An object of class "phylo". The output of 
        \code{parse.input.tree}, which preprocesses the
        observed tree.
    }
    \item{n.tips}{
        The number of tips in obs.tree.
    }
    \item{tip.heights}{
        A numeric vector containing the tip heights.
    }
    \item{tip.labels}{
        A character vector containing the tip labels as determined
        by the regular expression, regex, in \code{parse.labels}.
    }
    \item{sim.trees}{***
        A list of (number of particles) lists. Each list contains
        an integer from 1 to (number of particles). ***
    }
    \item{config}{
        The "smc.config" object given by the argument config.
    }
    \item{particles}{
        A (number of particles) x (number of parameters) matrix
        where each particle is a vector of model parameters.
    }
    \item{weights}{
        A vector of length (number of particles) containing NA
        to hold weights of particles.
    }
    \item{new.weights}{
        A vector of length (number of particles) containing NA
        to hold new weights of particles.
    }
    \item{dists}{
        A (number of samples) x (number of particles) matrix
        containing NA to hold distances from kernel.
    }
    \item{epsilon}{
        A numeric vector containing the current tolerance. The
        largest normalized floating-point number of the machine
        R is running on.
    }
    \item{accept}{
        A scalar vector containing the number of accepted proposals.
    }
    \item{alive}{
        A scalar vector containing the number of live particles.
    }
    \item{seed}{
        Contains NA if no seed is set.
    }
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \blod{22}(5), 1009--1020.
}

\author{
    Art Poon and Mathias Renaud
}

\note{
%%  ~~further notes~~
}

\seealso{
    \code{\link{parse.input.tree}} for the processing of obs.tree.
    \code{\link{smc.config}} for the assembly of the smc.config object.
    \code{\link{parse.labels}} for the processing of tree tip labels.
    \code{\link{print.smc.workspace}} for printing the result of \code{init.workspace}.
    \code{\link{run.smc}} for the use of \code{init.workspace}.
}

\examples{
## Set observed tree:
t1 <- read.tree(text="(A:0.1,B:0.2):0;")

## Set SMC configurations:
config <- list(
     params=c('N'),
     nparticle=10,
     nsample=3,
     decay.factor=1.0,
     rbf.variance=2.0,
     sst.control=1,
     norm.mode='NONE'
)
class(config) <- 'smc.config'

ws <- init.workspace(t1, config, regex=NA)
}

\keyword{SMC}
\keyword{Workspace}
\keyword{smc.workspace}