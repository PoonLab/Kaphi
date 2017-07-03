\name{run.smc}
\alias{run.smc}

\title{
    Run Sequential Monte Carlo
}

\description{
    \code{run.smc()} performs a Sequential Monte Carlo (SMC) for
    Approximate Bayesian Computation (ABC).
}

\usage{
    run.smc(ws, trace.file = "", regex = NA, seed = NA, nthreads = 1, verbose = FALSE, ...)
}

\arguments{
    \item{ws}{
        An object of class 'smc.workspace'.
    }
    \item{trace.file}{
        (optional) path to a file to write outputs. The default
        value of '' prints the output to the screen.
    }
    \item{regex}{
        (optional) A regular expression to pass to \code{parse.labels()}
        via \code{init.workspace()}. Default is NA.
    }
    \item{seed}{
        (optional) integer to set random seed, passed to \code{simulate.trees()}.
        Default is NA.
    }
    \item{nthreads}{
        (optional) for running on multiple cores. Default is 1.
    }
    \item{verbose}{
        If TRUE, additional information is printed to the screen
        with each iteration. Default is FALSE.
    }
    \item{\dots}{
        Additional arguments to pass to config@generator via \code{simulate.trees()}.
    }
}

\details{}

\value{
    A list of two lists: the final workspace and the results from the run.
    \item{workspace}{
      An smc.workspace object.
    }
    \item{result}{
      A list of returned values:
      \itemize{
        \item{niter: }{number of iterations.
        }
        \item{theta: }{
            a list of matricies containing the theta parameter values
            for each particle at each iteration.
        }
        \item{weights: }{
            a list of vectors containg the weight of each particle at
            each iteration.
        }
        \item{accept.rate: }{
            a vector containing the acceptance rate at each iteration.
        }
        \item{epsilons: }{
            a vector containing the value of epsilon at each iteration.
        }
      }  
    }
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}(5), 1009--1020.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{initialize.smc}}
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

## Initialize SMC
ws <- initialize.smc(ws)

result <- run.smc(ws)
}

\keyword{SMC}
\keyword{ ~kwd2 }