\name{plot.smc.config}
\alias{plot.smc.config}
\title{
Display Prior Distributions
}
\description{
Plots the different prior distributions within the Kaphi sequential Monte Carlo configuration settings. Multiple plots are generated, one plot each for every sample from the prior distribution.
}
\usage{
plot.smc.config(config, nreps=1000, numr=1, numc=1)
}
\arguments{
  \item{config}{
    A list of settings configured based on YAML input from the user. A \code{smc.config} object.
}
  \item{nreps}{
    Number of replicates. Default is 1000. 
  }
  \item{numr}{
    Number of rows of prior distributions to be displayed at one time. Default is 1.
  }
  \item{numc}{
    Number of columns of prior distributions to be displayed at one time. Default is 1.
  }
}
\note{
User will be prompted to advance onto the next plot in the graphical window by using the <Return> key.
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
