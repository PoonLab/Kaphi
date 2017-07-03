\name{prior.density}
\alias{prior.density}
\title{
Retrieve Prior Density
}
\description{
Retrieves density of the prior distribution. 
}
\usage{
prior.density(config, theta)
}
\arguments{
  \item{config}{
    A list of settings configured based on YAML input from the user. A \code{smc.config} object.
}
  \item{theta}{
    A named vector of the configuration parameters.
  }
}
\value{
Returns density of the prior distribution.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{sample.priors}}, \code{\link{propose}}, \code{\link{proposal.density}}
}
\examples{
require(Kaphi, quietly=TRUE)

#load configuration file
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-coalescent.yaml')    # Ne.tau parameter has prior distribution of: rlnorm(n=1,meanlog=5,sdlog=2)

theta <- c(Ne.tau=1000)                 # say this is the true value
prior.density(config, theta)
}