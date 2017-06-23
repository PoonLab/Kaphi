\name{tree.stats}
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
\alias{avg.unbalance}
\alias{internal.terminal.ratio}
\title{Compute Tree Statistics}
\description{
\code{sackin} is a wrapper function that returns the sum of node depths (branch lengths) from each tip to the root. \cr\cr
\code{colless} is a wrapper function that returns an index between 0 and 1 to represent the level of imbalance in the specified tree.\cr\cr
\code{cophenetic} is a wrapper function that returns the sum of most-recent common ancestor depths for each pair of tips.\cr\cr
\code{pybus.gamma} is a wrapper function that calculates and returns Pybus' gamma-statistic.\cr\cr
\code{ladder.length} is a wrapper function that computes and returns the maximum ladder length of the specified tree. \cr\cr
\code{IL.nodes} is a wrapper function that returns the number of internal nodes which have exactly one leaf child. \cr\cr
\code{tree.width} is a wrapper function that returns the maximum width of each tree divided by its maximum depth.\cr\cr
\code{max.delta.width} is a wrapper function that returns the absolute value of the maximum difference in widths (number of nodes) between two levels in a tree.\cr\cr
\code{n.cherries} is a wrapper function that returns the number of internal nodes with two child descendants as leaves in a tree.\cr\cr
\code{prop.unbalanced} is a wrapper function that finds and returns the proportion of unbalanced subtrees in the specified tree.\cr\cr
\code{avg.unbalance} is a wrapper function that computes and returns the average ratio of unbalanced subtrees in the specified tree.\cr\cr
\code{internal.terminal.ratio} is a wrapper function that computes and returns the ratio of internal to terminal branch lengths of the tree.
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
avg.unbalance(x)
internal.terminal.ratio(x)
}
\arguments{
  \item{x}{
    An object of class "phylo". Tree must be ultrameric in the case of Pybus' gamma-statistic.
  }
  \item{use.branch.lengths}
    If 0, all branches are treated as if they had unit length.
}
\details{
All wrapper functions call the associated function in the C extension \code{src/kernel.c}  \cr\cr
Sackin's imbalance is a summary metric used to summarize the topology of the tree. It is measured as the average length of the paths from the leaves to the root of the tree.\cr\cr  
The Colless imbalance is a normalized measure of the asymmetry of a rooted full binary tree. A completely asymmetric tree will have an imbalance of 1 and a symmetric tree will have an imbalance of 0.\cr\cr
The \emph{total cophenetic index} is the sum over all pairs of different leaves at the depth of their least common ancestor.  \cr\cr
Pybus' gamma-statistic summarizes the information contained in the inter-node intervals of a phylogeny. The tree should be ultrameric, but this isn't checked by the function, so if it is not ultrameric, the returned result may not be meaningful. \cr\cr
Ladder length is defined as the maximum path length in the tree, in untis of number of branches (edges, not the branch length). In other terms, it is the max number of connected internal nodes with a single leaf descendant, divided by the number of leaves in the tree. Also known as the trunk length.\cr\cr
The "IL" portion of the nodes is the proportion of internal nodes with a single leaf descendant.  \cr\cr
The width of a tree is described as the maximum number of nodes in the tree in terms of the depth of the tree. The "depth" of a node in a tree is the number of edges between itself and the tree's root, also known as the level of the tree. The "width" of the tree at a depth d is defined as the number of nodes in the tree with depth d. 
The maximum difference in widths of the tree is another statistic that reflects the absolute value of the maximum difference in number of nodes between two levels in a tree. \cr\cr
A cherry configuration is an internal node whose children are both leaves.\cr\cr
An unbalanced tree would have different numbers of descending tips on the left and right sides of the internal node.  \cr\cr
The average ratio of unbalanced substrees in the specified tree is defined as the average of min(L, R) / max(L, R) where L (with respect to R) is the number of tips descending from the left subtree (with respect to the right subtree).  \cr\cr
The ratio of itnernal to terminal branch lengths is described as the ratio of the average internal branch length to the average terminal branch length.
}
\references{
Colijn, C., & Gardy, J. (2014). Phylogenetic tree shapes resolve disease transmission patterns. \emph{Evolution, Medicine, & Public Health}, \bold{1}, 96--108.\cr\cr
Mir, A., Rossello F., & Rotger, L. (2013). A new balance index for phylogenetic trees. \emph{Mathematical Biosciences}, \bold{241.1}, 125--136.\cr\cr
Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proceedings of the Royal Society of London. Series B. Biological Sciences}, \bold{267}, 2267--2272.\cr\cr
Norstrom, M. M., Prosperi, M. C., Gray, R. R., et al. (2012). Phylotempo: a set of R scripts for assessing and visualizing temporal clustering in genealogies inferred from serially sampled viral sequences. \emph{Evolutionary Bioinformatics; London}, \bold{8}, 261.
}
\author{
Art Poon and Tammy Ng
}
\seealso{
\code{\link{nLTT}}
}
\examples{
require(Kaphi, quietly=TRUE)

#create tree in Newick format
tree <- read.tree(text="(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0;")

sackin(tree, TRUE)   # using branch lengths = TRUE
colless(tree)
cophenetic(tree)
pybus.gamma(tree)
ladder.length(tree)
tree.width(tree)
IL.nodes(tree)
max.delta.width(tree)
n.cherries(tree)
prop.unbalanced(tree)
avg.unbalance(tree)
internal.terminal.ratio(tree)
}