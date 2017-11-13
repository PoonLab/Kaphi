\name{MAST}
\alias{MAST}

\title{
Maximum Agreement SubTree distance metric
}
\description{
Topology-only distance that measures the similarity between two trees as the number of tips in the largest subtree identical between the two trees. Converted to a distance by subtracting the similarity score from the tree size.
}
\usage{
MAST(x, y)
}

\arguments{
  \item{x}{
    A bifurcating tree object of class 'phylo'
  }
  \item{y}{
    A bifurcating tree object of class 'phylo' of the same size (has the same number of tips) as tree \code{x}
  }
}

\value{
MAST score of type \code{numeric} between two given trees.
}
\references{
Gordon, A.D., & Tomassine, R. (1980). On the assessment and comparison of classifications, \emph{Analyse de donnes et informatique}, pg 149--160. \cr\cr
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214.
}
\author{
Tammy Ng
}

\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

MAST(t6,t7)
}