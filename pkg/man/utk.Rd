\name{utk}
\alias{utk}
\title{
Unlabelled tree-shaped kernel
}
\description{
Convenience wrapper for an unlabelled tree shape kernel.
}
\usage{
utk(x, y, config)
}
\arguments{
  \item{x}{
    An object of class "phylo"
    }
  \item{y}{
    An object of class "phylo"
  }
  \item{config}{
    A cached object of self-kernel
  }
}
\value{
Returns a tree kernel with labelled tree tips.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{tree.kernel}}
}
\examples{
require(Kaphi, quietly=TRUE)

#load config
setwd('~/git/Kaphi')
config <- load.config('pkg/examples/example-compartmental.yaml')

#create tree in Newick format
t1 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")
t2 <- read.tree(text="(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0;")
utk(t1, t2, config)
}
