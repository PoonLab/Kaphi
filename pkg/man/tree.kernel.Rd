\name{tree.kernel}
\alias{tree.kernel}
\alias{utk}
\title{
Tree shape kernel function
}
\description{
A kernel function maps objects into a high dimensional feature space where the inner product between objects can be computed efficiently.     
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
    Decay factor is used in the kernel function to penalize larger subset tree matches. Without penalities, large phylogenies will have enormous values along the diagonal of the kernel matrix because it will always match itself exactly, so then it washes out any other patterns when comparing with other phylogenies. 
  }
  \item{sigma}{
    Variance of Gaussian radial basis function of branch lengths.
  }
  \item{rho}{
    Subset tree control parameter. An integer value of 0 to represent a subtree kernel, and an integer value of 1 to represent a pure subset tree kernel. Default is 1.0.
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
  Moschitti, A. (2006). Making tree kernels practical for natural language learning. \emph{Association for Computational Linguistics - European Chapter}. \bold{24}, 113--120.
}
\author{
%%  ~~who you are~~
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
