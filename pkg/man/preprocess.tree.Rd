\name{preprocess.tree}
\alias{preprocess.tree}
\alias{PREPROCESS.TREE}
\alias{rescale.tree}
\alias{tree.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocess a Tree
}
\description{
Processes the tree to ensure that the class of the tree is of type character or an phylo object (package=ape) before reorganizing the internal structure of the tree to get the ladderized effect when plotted. Rescales the tree according to specified configuration's normal mode, and caches a self-kernel score. 
}
\usage{
  preprocess.tree(x, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class phylo
  }
  \item{config}{
    An cached object of the self-kernel?
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the tree that has been ladderized and rescaled, with a self-kernel score calculated and cached into the tree.
}
\author{
Author(s) go here
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tree.kernel}}
\code{\link{rescale.tree}}
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
