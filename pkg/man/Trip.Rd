\name{Trip}
\alias{Trip}

\title{
Triplet distance metric
}
\description{
Topology-only distance that enumerates all possible triplets of the tips of two given rooted trees. Determines topology of that triplet in A and B and scores 1 point for each triplet whose topology differs, and sums this score over the triplets. This implementation divides the result by the total number of triplets.
}
\usage{
Trip(x, y)
}

\arguments{
  \item{x}{
    A rooted tree object of class 'phylo'
  }
  \item{y}{
    A rooted tree object of class 'phylo'
  }
}

\value{
Triplets score of type \code{numeric} between two given trees.
}
\references{
Critchlow, D.E., Pearl, D.K., & Qian, C. (1996). The triples distance for rooted bifurcating phylogenetic trees. \emph{Systematic Biology}, \bold{45}, 323--334. \cr\cr
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214.
}
\author{
Mathias Renaud and Tammy Ng
}

\seealso{
  \code{\link{TripL}}
}
\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

Trip(t6,t7)
}

