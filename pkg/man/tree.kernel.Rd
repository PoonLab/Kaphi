\name{tree.kernel}
\alias{tree.kernel}
\alias{rescale.tree}
\alias{preprocess.tree}
\alias{to.newick}
\alias{utk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tree shape kernel functions
}
\description{
A kernel function maps objects into a high dimensional feature space where the inner product between objects can be computed efficiently.  
Kaphi uses a   
\code{tree.kernel} is a wrapper function that calls the associated function in the C extension \code{src/kernel.c} to compute a tree shape kernel function. 
}
\usage{
tree.kernel(x, y, lambda, sigma, rho=1.0, normalize=0, label1=NA, label2=NA, gamma=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class "phylo"
}
  \item{y}{
    An object of class "phylo"
  }
  \item{lambda}{
    Decay factor of type ...
  }
  \item{sigma}{
    (define) RBF variance parameter of type ...
  }
  \item{rho}{
    (define) SST control parameter. An integer value of 0 to represent a subtree kernel, and an integer value of 1 to represent a subset tree kernel. Default is 1.0.
  }
  \item{normalize}{
    Normalizes kernel score by sqrt(k(x,x) * k(y,y)). Default is 0 (but what does 0 mean?)
  }
  \item{label1}{
    Argument for a label to be assigned to the tips of tree x, an object of class "phylo". Default is NA.
  }
  \item{label2}{
    Argument for a label to be assigned to the tips of treey y, an object of class "phylo". Default is NA.
  }
  \item{gamma}{
    (Not sure what this parameter does)
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns results of the associated C extension function that computes a tree shaped kernel function
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
\code{\link{to.newick}}
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
