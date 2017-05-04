\name{get.times}
\alias{get.times}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Generate vector of time points}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.times(t0, t1, time.pts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{t1}{
%%     ~~Describe \code{t1} here~~
}
  \item{time.pts}{
%%     ~~Describe \code{time.pts} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (t0, t1, time.pts) 
{
    times1 <- seq(t0, t1, length.out = time.pts)
    tdelta <- times1[2] - times1[1]
    if ((t0 - tdelta) <= t0) {
        times0 <- c()
    }
    else {
        times0 <- seq(t0, t0 - tdelta, by = tdelta)
    }
    times <- unique(c(times0, times1))
    return(times)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
