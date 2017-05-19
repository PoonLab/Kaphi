\name{prior.density}
\alias{prior.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve Prior Density
}
\description{
Retrieves density of the prior distribution. 
}
\usage{
prior.density(config, theta)
}
%- maybe also 'usage' for other objects documented here.
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
%%  ~~who you are~~
}
\seealso{
\code{\link{sample.priors}}, \code{\link{propose}}, \code{\link{proposal.density}}
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
