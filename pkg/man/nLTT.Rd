\name{nLTT}
\alias{nLTT}
\title{
  Compute the normalized lineages-through-time (nLTT) statistic
}
\description{
  Slightly modified from the version in \cite{Janzen et al.} using the trapezoid rule instead of the rectangle method. This is a wrapper function that calls the associated function in the C extension \code{src/kernel.c} to make available other tree shape similarity measures and statistics.
}
\usage{
nLTT(x, y)
}
\arguments{
  \item{x}{
    An object of class "phylo". The first tree to compare.
}
  \item{y}{
    An object of class "phylo". The second tree to compare.
  }
}
\value{
Returns the nLTT statistic.
}
\references{
Janzen, T., Hohna, S., and Etienne, R. S. (2015). Approximate Bayesian Computation of diversification rates from molecular phylogenies: introducing a new efficient summary statistic, the nLTT. \emph{Methods in Ecology and Evolution}, \bold{6:5}, 566--575.
}
\author{
Art Poon and Tammy Ng
}
\details{
The normalized lineages-through-time statistic is primarily suited for comparing trees rather than providing information about a single tree. May be used as a substitute for the likelihood. The LTT curve shows the number of lineages in the phylogeny over time. Normalizing two phylogenies in terms of number of lineages and time provides a level of comparison. 
}

\seealso{
\code{\link{tree.stats}}
}
\examples{
####issue #49 not resolved yet
}

