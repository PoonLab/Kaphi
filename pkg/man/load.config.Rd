\name{load.config}
\alias{load.config}
\title{
Loads Configuration of..?
}
\description{
Loads file and sets up configuration with the default settings, including sequential Monte-Carlo (SMC) settings and kernel settings. 
}
\usage{
load.config(filename)
}
\arguments{
  \item{filename}{
  File must be a YAML string to be loaded into the configuration.
  }
}
\details{
YAML string will return a sequence that is not uniform, therefore a list is returned. 
}
\value{
A list is returned of items stored in \code{config}, an \code{smc.config} object.
  \item{params}{
    Contains the names of the priors.
}
  \item{priors}{
    Name(s) of the prior distribution(s)
  }
  \item{prior.densities}{
    The prior probability distribution in the form of a \code{list}.
  }
  \item{proposals}{
    Name(s) of the proposal distribution(s). Variable set in \code{proposals} must be a \code{list} and a subset of \code{priors}
  }
  \item{proposal.densities}{
    The proposal distribution that is proportional to the target distribution. Generated from making the sequence of samples into a random walk in the form of a \code{list}. The proposal distribution \strong{g(x'|x)} is the conditional probability of proposing a state \strong{x'} given \strong{x}.
  }
  \item{model}{
    Model used for generating proposal distribution from given prior distribution?
  }
  \item{nparticle}{
    Number of particles, or Monte-Carlo samples. Default is 1000.
  }
  \item{nsample}{
    Number of samples per particle. Default is 5. 
  }
  \item{ess.tolerance}{
    Explained sum of squares (ESS) is a quantity describing how well the model represents the data being modelled. It measures how much variation there is in the modelled values compared to the total sum of squares (variation in the observed data) as well as the residual sum of squares (variation in modelling errors). Default is 1.5.
  }
  \item{final.epsilon}{
    Final epsilon value. Default is 0.01.
  }
  \item{final.accept.rate}{
    Final acceptance rate. Default is 0.015.
  }
  \item{quality}{
    Quality of the distribution? Default is 0.95.
  }
  \item{step.tolerance}{
    The tolerance level decreases by this amount at subsequent algorithmic time steps until it reaches \code{final.epsilon}. If \code{step.tolerance} decreases too slowly then the algorithm will be too computationally intensive. If tolerance levels decrease too fast, the algorithm can perform poorly. Default is 1e-5.
  }
  \item{decay.factor}{
    Decay factor is used in the kernel function to penalize larger subset tree matches. Default is 0.2.
  }
  \item{rbf.variance}{
    Variance of Gaussian radial basis function of branch lengths. Default is 2.0.
  }
  \item{sst.control}{
    Subset tree control parameter. An integer value of 0 represents a subtree kernel, and an integer value of 1 to represent a pure subset tree kernel. Default is 1.0.
  }
  \item{norm.mode}{
    Normalization mode. Default is 'MEAN'.
  }
}
\references{
Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}, 1009--1020.
}
\author{
%%  ~~who you are~~
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
