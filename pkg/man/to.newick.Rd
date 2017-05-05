\name{to.newick}
\alias{to.newick}
\alias{write.tree}
\alias{read.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks and Converts Tree to Newick Format
}
\description{
Checks that the tree argument is an phylo object (package = ape) or a character string in standard Newick (New Hampshire) format before writing the tree file into parenthetic format using Newick format. 
}
\usage{
to.newick(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class "phylo" or a character string in standard Newick format.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector of mode character if file = "".
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Author(s)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{write.tree}}
\code{\link{read.tree}}
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
