\name{sackin}
\alias{sackin}
\alias{colless}
\alias{cophenetic}
\alias{pybus.gamma}
\alias{ladder.length}
\alias{IL.nodes}
\alias{tree.width}
\alias{max.delta.width}
\alias{n.cherries}
\alias{prop.unbalanced}
\alias{avg.unbalanced}
\alias{internal.terminal.ratio}
\title{
Compute Tree Statistics
}
\usage{
sackin(x, use.branch.lengths=FALSE)
colless(x)
cophenetic(x. use.branch.lengths=FALSE)
pybus.gamma(x)
ladder.length(x)
IL.nodes(x)
tree.width(x)
max.delta.width(x)
n.cherries(x)
prop.unbalanced(x)
avg.unbalanced(x)
internal.terminal.ratio(x)
}
\arguments{
  \item{x}{
    An object of class "phylo". Tree must be ultrameric in the case of Pybus' gamma statistic.
  }
  \item{use.branch.lengths}
    If 0, all branches are treated as if they had unit length.
}

\subtitle{Sackin's Index}
\description{
Sackin's imbalance is a summary metric used to summarize the topology of the tree. It is measured as the average length of the paths from the leaves to the root of the tree. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} returns the sum of node depths (branch lengths) from each tip to the root.
}

\subtitle{Colless' Index}
\description{
The Colless imbalance is a normalized measure of the asymmetry of a rooted full binary tree. A completely asymmetric tree will have an imbalance of 1 and a symmetric tree will have an imbalance of 0. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} and returns Colless' index.
}

\subtitle{Cophenetic Index}
\description{
The total cophenetic index is the sum over all paris of different leaves of the depth of their least common ancestor. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} and returns the sum of most-recent common ancestor depths for each pair of tips.
\references{
Mir, A., Rossello F., & Rotger, L. (2013). A new balance index for phylogenetic trees. \emph{Mathematical Biosciences}, \bold{241.1}, 125--136.
}

\subtitle{Pybus' gamma-statistic}
\description{
Gamma-statistic summarizes the information contained int he inter0node intervals of a phylogeny. The tree should be ultrametric, but this isn't checked by the function, so if it is not, the returned result may not be meaningful. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} and returns the gamma statistic, a numeric vector of length one.
}
\references{
Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proceedings of the Royal Society of London. Series B. Biological Sciences}, \bold{267}, 2267--2272.
}

\subtitle{ladder.length}
\description{
Computes the maximum ladder length. This is the maximum path length in the tree, in units of number of branches (not the branch length). In other terms, it is the max number of connected internal nodes with a single leaf descendant, divided by the number of leaves in the tree. Also known as the trunk length. A wrapper function that calls the associated function in the C extension \code{src/kernel.c}.
}
\references{
Colijn, C., & Gardy, J. (2014). Phylogenetic tree shapes resolve disease transmission patterns. \emph{Evolution, Medicine, & Public Health}, \bold{1}, 96--108.
}

\subtitle{IL.nodes}
\description{
The proportion of internal nodes with a single leaf descendant. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} that returns the number of internal nodes which have exactly one leaf child. 
}

\subtitle{tree.width}
\description{
Computes the maximum number of nodes at any depth (not considering branch lengths). The 'depth' of a node in a tree is the number of edges between itself and the tree's root. The 'width' of a tree at a depth d is defined as the number of nodes with depth d. A wrapper function that calls the associated function in the C extension \code{src/kernel.c} that returns the maximum width of each tree divided by its maximum depth.
}

\subtitle{max.delta.width}
\description{
Reflects the maximum absolute difference in widths from one depth to the next, over all depths d in the tree.
}

\value{

}
\references{

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

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
