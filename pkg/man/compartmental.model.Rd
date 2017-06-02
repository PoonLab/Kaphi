\name{compartmental.model}
\alias{compartmental.model}
\title{
Simulate Trees with specified Compartmental Model
}
\description{
Calculates a numerical solution of an ODE system specified by a given compartmental model and returns an ape multiPhylo object of simulated trees.
}
\usage{
compartmental.model(theta, nsim, tips, model='sir.nondynamic', seed=NA, fgyResolution=500, integrationMethod='adams')
}

\arguments{
  \item{theta}{
    Vector containing parameter values for the given model. Should contain
  }
  \item{nsim}{
    Number of replicate trees to simulate
  }
  \item{tips}{
    Number of tips of zero height (\code{integer}) OR vector of tip heights (\code{vector})
  }
  \item{model}{
    Specified compartmental model. \code{'sir.nondynamic'} is a susceptible-infectious-recovered model without vital dynamics (births and deaths). \code{'sir.dynamic'} is a susceptible-infectious-recovered model with vital dynamics and constant population. \code{'sis'} is a susceptible-infectious-susceptible model with vital dynamics. \code{'seir'} is a susceptible-exposed-infectious-recovered model with vital dynamics where birth rate is equal to death rate. Default is \code{'sir.nondynamic'}.
  }
  \item{seed}{
    Setting seed for pseudorandom generator. Default is NA.
  }
  \item{fgyResolution}{
    Time resolution of ODE solution
  }
  \item{integrationMethod}{
    The integrator to use. Default integrator used is \code{adams}.
  }
}
\details{
The integration method used is either a function that performs integration, or a list of class rkMethod, or a string ("lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" ,"iteration"). Options "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" are the backward differentiation formula, the BDF with diagonal representation of the Jacobian, the (explicit) Adams and the implicit Adams method, and the implicit Adams method with diagonal representation of the Jacobian respectively (see details). 
}
\value{
An ape \code{'multiPhylo'} object of simulated trees.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{speciation.model}} , \code{\link{const.coalescent}}
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

