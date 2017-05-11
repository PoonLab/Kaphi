\name{load.config}
\alias{load.config}
\title{
Loads Configuration of..?
}
\description{
Loads configuration with the default settings, including sequential Monte-Carlo (SMC) settings and kernel settings. 
}
\usage{
load.config(params=NA,
            priors=list(),
            prior.densities=list(),
            proposals=list(),
            proposal.densities=list(),
            model=NA,

            # SMC settings
            nparticle=1000,
            nsample=5,
            ess.tolerance=1.5,
            final.epsilon=0.01,
            final.accept.rate=0.015,
            quality=0.95,
            step.tolerance=1e-5,

            # kernel settings
            decay.factor=0.2,
            rbf.variance=2.0,
            sst.control=1.0,
            norm.mode='MEAN')
}
\arguments{
  \item{params}{
%%     ~~Describe \code{x} here~~
}
  \items{priors}{
  
  }
  \items{prior.densities}{
  
  }
  \items{proposals}{
    Variable set in \code{proposals} must be a \code{list} and a subset of \code{priors}
  }
  \items{proposal.densities}{
  
  }
  \items{model}{
    
  }
  \items{nparticle}{
    Number of particles, or Monte-Carlo samples. Default is 1000.
  }
  \items{nsample}{
    Number
  }
  \items{ess.tolerance}{
  
  }
  \items{final.epsilon}{
  
  }
  \items{final.accept.rate}{
  
  }
  \items{quality}{
  
  }
  \items{step.tolerance}{
    The tolerance level decreases by this amount at subsequent algorithmic time steps until it reaches \code{final.epsilon}. 
  }
  \items{decay.factor}{
    Decay factor is used in the kernel function to penalize larger subset tree matches. Default is 0.2.
  }
  \items{rbf.variance}{
    Variance of Gaussian radial basis function of branch lengths. Default is 2.0.
  }
  \items{sst.control}{
    Subset tree control parameter. An integer value of 0 represents a subtree kernel, and an integer value of 1 to represent a pure subset tree kernel. Default is 1.0.
  }
  \items{norm.mode}{
    Default is 'MEAN'.
  }
}
\details{
If \code{step.tolerance} decreases too slowly then the algorithm will be too computationally intensive. If tolerance levels decrease too fast, the algorithm can perform poorly.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}, 1009--1020.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
