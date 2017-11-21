\name{Node.dist}
\alias{Node.dist}

\title{
Nodal/Path distance metric
}
\description{
Topology-only distance metric that considers each possible pair of typs in turn and counts the number of nodes traversed in the minimal path from tip A to tip B in a tree. The distance between tree \code{x} and tree \code{y} is the sum of the differences in minimal path lengths.
}
\usage{
Node.dist(x, y, k=0.65)
}

\arguments{
  \item{x}{
    A tree object of class 'phylo'
  }
  \item{y}{
    A tree object of class 'phylo' of the same size (has the same number of tips) as tree \code{x}
  }
  \item{k}{
    Power that the differences in minimal path lengths will be raised to before summing. Default is 0.65.
  }
}

\value{
Nodal/path distance score of type \code{numeric} between two given trees. Sum of the differences in minimal path lengths between tree \code{x} and tree \code{y}
}
\references{
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214. \cr\cr
Williams, W.T., & Clifford, H.T. (1971). On the comparison of two classifications of the same set of elements. \emph{Taxon}, \bold{20}, 519--522.
}
\author{
Tammy Ng
}

\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

Node.dist(t6,t7, k=1)
}