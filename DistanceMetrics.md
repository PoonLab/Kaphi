#### Valid Distance Metrics
| Tree Statistic | Package | Description |
|----------------|---------|-------------|
| kernel.dist*| Kaphi |The kernel distance between two trees.|
| ~~nLTT~~| ~~Kaphi~~|*Does not currently work*|
| sackin | Kaphi |Sackin index.|
| colless | Kaphi |Colless imbalance number.|
| cophenetic.index | Kaphi |Cophenetic index.|
| ladder.length | Kaphi |Max ladder length.|
| IL.nodes | Kaphi |Number of internal nodes with one leaf.|
| tree.width | Kaphi |Max width divided by max depth.|
| max.delta.width | Kaphi |Max difference in in width between two levels.|
| n.cherries | Kaphi |Number of cherries(node with two leaves).|
| prop.unbalanced | Kaphi |Proportion of unbalanced subtrees.|
| avg.unbalance | Kaphi |Average ratio of unbalanced subtrees.|
| pybus.gamma | Kaphi |Pybus' gamma statistic|
| internal.terminal.ratio | Kaphi |Ratio of internal to terminal branches.|
| cophenetic.phylo.met | Kaphi |Wrapper of `ape::cophenetic.phylo`, computes pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.|
| dist.nodes.met | Kaphi |Wrapper of `ape::dist.nodes`, similar to `cophenetic.phylo`, but includes internal nodes.|
| dist.topo | ape |Topological distance between two trees using the method from Penny & Hendy (1985).|
| avgladder | phyloTop |Mean size of ladders in the tree.|
| getDepths.met | Kaphi |Wrapper of `phyloTop::getDepths`, originally returns a list of two vectors: `tipDepths` and `nodeDepths`, but given metric will output non-scalar values of type 'tips' or type 'nodes' specified by user.|
| pitchforks | phyloTop |Number of clades with three tips. |
| RF.dist | phangorn | Robinson-Foulds distance: number of clades not shared. |
| KF.dist | phangorn | RF with branch lengths: sums of differences in branch lengths. |
| path.dist | phangorn | Path distance metric: nodal distance ("NODE") with k=2. |
| Trip | Critchlow et al. (1996)| Proportion of triplets not shared between trees. |
| TripL | Kuhner & Yamato (2014)| Trip distance using branch lengths |
| MAST  | Gordon (1980) | Number of tips in the Maximum Agreement SubTree of given subtrees.  | 
| Align  | Nye et al. (2006)  | Calculates dissimilarity of of all one-to-one mapping of branches between 2 trees and determines most optimal mapping of branches  | 
| Sim  | Hein et al. (2005)  | Similiarity measure based on the probability that a point chosen randomly in A will be on a branch leading to the same set of tips as a point randomly chosen in B  | 
| Node  | Williams & Clifford (1971)  | Number of nodes traversed in mainimal path from one tip to another in trees A and B  | 
