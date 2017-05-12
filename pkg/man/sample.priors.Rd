\name{sample.priors}
\alias{sample.priors}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a named vector of sampled parameters from the prior distributions.
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
\code{\link{prior.density}}, \code{\link{propose}}, \code{\link{proposal.density}}
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
