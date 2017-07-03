\name{sample.priors}
\alias{sample.priors}
\title{
Sample the parameters from the prior distribution
}
\description{
Calls a function that samples parameters from the prior distributions. 
}
\usage{
sample.priors(config)
}
\arguments{
  \item{config}{
    List of settings configured based on YAML input from the user. A \code{smc.config} object.
}
}
\value{
Returns a named vector of sampled parameters from the prior distributions.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{prior.density}}, \code{\link{propose}}, \code{\link{proposal.density}}
}
\examples{
require(Kaphi, quietly=TRUE)

#load configuration file
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-coalescent.yaml')    # Ne.tau parameter has prior distribution of: rlnorm(n=1,meanlog=5,sdlog=2)

# sample priors of parameter(s)
theta <- sample.priors(config)
theta
}