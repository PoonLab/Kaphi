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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{prior.density}}, \code{\link{propose}}, \code{\link{sample.priors}}
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
