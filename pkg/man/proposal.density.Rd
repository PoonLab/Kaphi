\name{proposal.density}
\alias{proposal.density}
\title{
Retrieve Proposal Density
}
\description{
Retrieves density of the proposal distribution.
}
\usage{
proposal.density(config, theta, new.theta)
}
\arguments{
  \item{config}{
    A list of settings configured based on YAML input from the user. A \code{smc.config} object.
}
  \item{theta}{
    A named vector of the configuration parameters.
  }
  \item{new.theta}{
    A named vector of the parameters used for proposal distribution. A sanity check to ensure \code{theta} and \code{new.theta} have matched named vectors.
  }
}
\value{
Returns the density of the proposal distribution.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{prior.density}}, \code{\link{propose}}, \code{\link{sample.priors}}
}
\examples{
require(Kaphi, quietly=TRUE)

#load configuration file
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-coalescent.yaml')    # Ne.tau parameter has proposal distribution of: rnorm(n=1,mean=0,sd=50)

theta <- c(Ne.tau=1000)                 # say this is the true value
new.theta <- propose(config, theta)     # sample draw of parameter(s) of proposal distribution

proposal.density(config, theta, new.theta)
}