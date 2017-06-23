\name{propose}
\alias{propose}
\title{
Draw new values and Update proposal density
}
\description{
Draw new values from the proposal density and apply proposal values to update the parameter vector \code{theta}.
}
\usage{
propose(config, theta)
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
Returns \code{theta}, an updated vector of the configuration parameters.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{sample.priors}}, \code{\link{prior.density}}, \code{\link{proposal.density}}
}
\examples{
require(Kaphi, quietly=TRUE)

#load configuration file
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-coalescent.yaml')    # Ne.tau parameter has proposal distribution of: rnorm(n=1,mean=0,sd=50)

theta <- c(Ne.tau=1000)    #say this is the true value

# sample draw of parameter(s) of proposal distribution
new.theta <- propose(config, theta)
new.theta
}