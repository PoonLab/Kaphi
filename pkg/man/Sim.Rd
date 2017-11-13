\name{Sim}
\alias{Sim}

\title{
Similarity distance measure
}
\description{
Topology and branch lengths distance metric derived from a measure of similarity based on probability that a point randomly chosen in A will be on a branch leading to the same set of tips as a point randomly chosen in B. Scores are normalized by dividing similarities of A to B by the probability that the two points chosen at random in A will yield the same set fo tips. Converted to a distance metric by averaging the A/B and B/A similarities and subtracted from 1. 
}
\usage{
Sim(x, y)
}

\arguments{
  \item{x}{
    A tree object of class 'phylo' with branch lengths
  }
  \item{y}{
    A tree object of class 'phylo' with branch lengths
  }
}
\details{
Branch lengths are necessary, but the metric is not sensitive to their magnitudes, only to their ratios in the two trees. Two trees that differ only by a scaling factor in their branch lengths will test as identical.
}
\value{
Similarity score of type \code{numeric} between two trees.
}
\references{
Hein, J., Scierup, M.H., & Wiuf, C. (2005). Gene genealogies, variation and evolution: a primer in coalescent theory. \emph{Oxford University Press} \cr\cr
Kuhner, M.K., & Yamato, J. (2014). Practical Performance of Tree Comparison Metrics. \emph{Systematic Biology}, \bold{64}, 205--214.
}
\author{
Tammy Ng
}
\section{Warning }{Sim cannot be computed when either tree has a total branch length of zero, as the normalization step would divide be zero.}
\examples{
require(Kaphi, quietly=TRUE)
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")

Sim(t6,t7)
}