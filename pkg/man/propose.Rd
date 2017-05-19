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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sample.priors}}, \code{\link{prior.density}}, \code{\link{proposal.density}}
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
