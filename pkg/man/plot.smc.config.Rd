\name{plot.smc.config}
\alias{plot.smc.config}
\title{
Display Prior Distributions
}
\description{
Plots the different prior distributions within the Kaphi sequential Monte Carlo configuration settings. Multiple plots are generated, one plot each for every sample from the prior distribution.
}
\usage{
plot.smc.config(config, nreps=1000)
}
\arguments{
  \item{config}{
    A list of settings configured based on YAML input from the user. A \code{smc.config} object.
}
}
  \item{nreps}{
    Number of replicates. Default is 1000. 
  }
}
\note{
User has to advance onto the next plot using the right arrow key in the graphical window.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\seealso{
\code{\link{print.smc.config}}, \code{\link{sample.priors}}
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
