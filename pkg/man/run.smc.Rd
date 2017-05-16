\name{run.smc}
\alias{run.smc}

\title{
    Run Sequential Monte Carlo
}

\description{

}

\usage{
run.smc(ws, trace.file = "", regex = NA, seed = NA, nthreads = 1, verbose = FALSE, ...)
}

\arguments{
    \item{ws}{
        smc.workspace object.
    }
    \item{trace.file}{
        (optional) path to a file to write outputs.
    }
    \item{regex}{
        
    }
    \item{seed}{
        (optional) integer to set random seed.
    }
    \item{nthreads}{
        (optional) for running on multiple cores.
    }
    \item{verbose}{
%%     ~~Describe \code{verbose} here~~
    }
    \item{\dots}{
        additional arguments to pass to config@generator via simulate.trees().
    }
}

\details{}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}(5), 1009--1020.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }