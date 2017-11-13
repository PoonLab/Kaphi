\name{TripL}
\alias{TripL}

\title{
Triplet with branch lengths distance metric
}
\description{
Topology and branch lengths distance metric, computing the sum of discrepant branch lengths among triplets between two clock-like rooted trees.
}
\usage{
TripL(x, y, k=2)
}
\arguments{
  \item{x}{
    A rooted tree object of class 'phylo'
  }
  \item{y}{
    A rooted tree object of class 'phylo'
  }
  \item{k}{
    Power that the differences are raised to before summing. Default is 2.
  }
}
\details{
The length from the tips to the first node ('first-node lengths') from each tree are compared and the absolute values of their differences are summed.  Length from the first nodes to the second nodes ('second-node lengths') are similarly compared, and the results are summed across all triplets. Differences are raised to the power \emph{k} before summing.
}
\value{
Triplets with lengths score of type \code{numeric} between two given trees.
}
\references{
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214.
}
\author{
Mathias Renaud and Tammy Ng
}

\seealso{
\code{\link{Trip}}
}
\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

TripL(t6,t7)
}