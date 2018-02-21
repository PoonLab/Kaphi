\name{tree.kernel}
\alias{tree.kernel}
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
  \item{regexPattern}{
    Character string containing a regular expression to be matched in tree x and tree y tip labels.
  }
  \item{regexReplacement}{
    A replacement for matched \code{pattern} in tree x and tree y tip labels. Coerced to character if possible. This can include backreferences "\1" to "\9" to parenthesized subexpressions of \code{pattern}. If NA, all elements in the result corresponding to matches will be set to NA.
  }
  \item{gamma}{
    (Not sure what this parameter does)
  }
}
\value{
  Returns results of the associated C extension function that computes a tree shaped kernel function
}
\seealso{\code{\link{utk}}}
\references{
  Moschitti, A. (2006). Making tree kernels practical for natural language learning. \emph{Association for Computational Linguistics - European Chapter}. \bold{24}, 113--120.
}
\author{
Art Poon and Tammy Ng
}
\examples{
require(Kaphi, quietly=TRUE)

#create two trees in Newick format
t1 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")
t2 <- read.tree(text="(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0;")
lambda <- 2 
sigma <- 0.1

tree.kernel(t1, t2, lambda, sigma)
}
