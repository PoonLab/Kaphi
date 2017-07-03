\name{print.smc.config}
\alias{print.smc.config}
\title{
Print Kaphi SMC configuration
}
\description{
Prints the settings for the Kaphi sequential Monte Carlo configuration.
}
\usage{
print.smc.config(config, ...)
}
\arguments{
  \item{config}{
    A list of settings configured based on YAML input from the user. A \code{smc.config} object.
}
  \item{...}{
    Additional options ...??
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{plot.smc.config}}
}
\examples{
require(Kaphi, quietly=TRUE)
 
#load configuration file
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-coalescent.yaml')

print.smc.config(config)
# or it will print when you just write config variable
# config
}