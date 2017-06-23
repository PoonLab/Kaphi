\name{set.model}
\alias{set.model}
\title{
Set or Modifies Model in Configuration
}
\description{
Sets or modifies model in list of configurated settings, returning the reconfigured settings. 
}
\usage{
set.model(config, generator)
}
\arguments{
  \item{config}{
    The configuration of type \code{list} holding the model setting to be modified. A \code{smc.config} object.
}
  \item{generator}{
    Can either be a function name or an object.
  }
}
\value{
Returns a list of the reconfigured settings. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{load.config}}, 
}
\examples{
require(Kaphi, quietly=TRUE)
setwd('~/git/Kaphi')

config <- load.config('pkg/examples/example-coalescent.yaml')
config <- set.model(config, const.coalescent) 
}