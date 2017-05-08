\name{nLTT}
\alias{nLTT}
\title{
  Compute the normalized lineages through time (nLTT) statistic
}
\description{
  Slightly modified from the version in \cite{jansen2015approximate} using the trapezoid rule instead of the rectangle method. This is a wrapper function that calls the associated function in the C extension \code{src/kernel.c} to make available other tree shape similarity measures and statistics.
}
\usage{
nLTT(x, y)
}
\arguments{
  \item{x}{
    An object of class "phylo". The first tree to compare.
}
  \item{x}{
    An object of class "phylo". The second tree to compare.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the nLTT statistic
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
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
