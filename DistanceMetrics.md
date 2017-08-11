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
| cophenetic.phylo | ape |Pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.|
| dist.nodes | ape |Like `cophenetic.phylo`, but includes internal nodes.|
| dist.topo | ape |Topological distance between two trees using the method from Penny & Hendy (1985).|
| avgladder | phyloTop |Mean size of ladders in the tree.|
| getDepths | phyloTop |Returns a list of two vectors: `tipDepths` and `nodeDepths`.|
| pitchforks | phyloTop |Number of clades with three tips. |
