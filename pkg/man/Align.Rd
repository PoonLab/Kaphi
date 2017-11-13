\name{Align}
\alias{Align}

\title{
Tree alignment distance metric
}
\description{
Topology-only distance metric that considers all the ways in which a one-to-one mapping of branches in the two trees could be made. Calculates a dissimilarity score for the clades separated by each branch for each mapping, and takes the dissimilarity score of the optimal mapping as the distance between the trees.
}
\usage{
Align(x, y)
}
\arguments{
  \item{x}{
    A tree object of class 'phylo'
  }
  \item{y}{
    A tree object of class 'phylo' of the same size (has the same number of tips) as tree \code{x}
  }
}
\details{
Can be bifurcating or multifurcating trees. 
}
\value{
Alignment score of type \code{numeric} between two given trees.
}
\references{
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214. \cr\cr
Nye, T.M.W., Lio, P., & Gilks, W.R. (2006). A novel algorithm and web-based tool for comparing two alternative phylogenetic trees. \emph{Bioinformatics}, \bold{22}, 117--119.
}
\author{
Tammy Ng
}

\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

Align(t6,t7)
}