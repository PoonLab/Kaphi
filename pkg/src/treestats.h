/** \file treestats.h
 * \brief Functions for calculating metrics and summary statistics of trees.
 */

#ifndef TREESTATS_H
#define TREESTATS_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>

/* Blum, Michael GB, Olivier François, and Svante Janson. "The mean,
 * variance and limiting distribution of two statistics sensitive to
 * phylogenetic tree balance." The Annals of Applied Probability
 * (2006): 2195-2214.
 */
#define HARMONIC(n) ( M_EULER + gsl_sf_psi((n) + 1) )
#define SACKIN_YULE(n) ( 2 * (n) * (HARMONIC(n) - 1) )
#define COLLESS_YULE(n)  ( (n) * log(n) + (n) * (M_EULER - 1 - log(2)) )
#define COPHENETIC_YULE(n)  ( (n) * ((n) - 1) - 2 * (n) * (HARMONIC(n) - 1) )

/** Calculate the tree kernel.
 *
 * Uses the fast algorithm from \cite moschitti2006making.
 * See \cite poon2015phylodynamic for the meanings of the parameters
 * decay_factor and rbf_variance.
 *
 * \param[in] t1,t2 trees to compare
 * \param[in] decay_factor decay factor in [0, 1] penalizing large matches
 * \param[in] rbf_variance variance of Gaussian radial basis function of branch lengths
 * \param[in] sst_control between 0 and 1, where 0 is a pure subtree kernel and
 * 1 is a pure subset tree kernel (see \cite moschitti2006making)
 */
double kernel(const igraph_t *t1, const igraph_t *t2, double decay_factor,
        double rbf_variance, double sst_control);

/** Compute the normalized lineages-through-time statistic.
 *
 * Slightly modified from the version in \cite janzen2015approximate, using the
 * trapezoid rule instead of the rectangle method.
 *
 * \param[in] t1 first tree to compare
 * \param[in] t2 second tree to compare
 * \return the nLTT statistic
 */
double nLTT(const igraph_t *t1, const igraph_t *t2);

/** Compute Sackin's index.
 *
 * \param[in] t tree to compute Sackin's index for
 * \param[in] use_branch_lengths if 0, treat all branches as if they had unit
 * length
 * \return the average path length from tips to the root
 */
double sackin(const igraph_t *t, int use_branch_lengths);

/** Compute Colless' index.
 *
 * \param[in] t tree t compute Colless' index for
 * \return Colless' index
 */
int colless(const igraph_t *t);

/** Compute the total cophenetic index.
 *
 * See Mir, Arnau, and Francesc Rosselló. "A new balance index for phylogenetic
 * trees." Mathematical biosciences 241.1 (2013): 125-136.
 *
 * \param[in] tree tree to compute index for
 * \param[in] use_branch_lengths if 0, treat all branches as if they had unit length
 * \return the sum of most-recent common ancestor depths for each pair of tips
 */
double cophenetic(const igraph_t *tree, int use_branch_lengths);

/** Compute the maximum ladder length.
 *
 * This is the maximum path length in the tree, in units of number of branches
 * (not branch length). Also known as the trunk length.  See \cite
 * colijn2014phylogenetic.
 *
 * \param[in] tree tree to compute ladder length for
 * \return the maximum number of branches from the root to a tip
 */
int ladder_length(const igraph_t *tree);

/** Compute the number of IL nodes.
 *
 * IL nodes are those with exactly one leaf child. See \cite
 * colijn2014phylogenetic.
 *
 * \param[in] tree the tree to compute the number of IL nodes for
 * \return the number of internal nodes which have exactly one leaf child
 */
int il_nodes(const igraph_t *tree);

/** Compute the width of a tree.
 *
 * The width is the maximum number of nodes at any depth (not
 * considering branch lengths). See \cite colijn2014phylogenetic.
 *
 * \param[in] tree tree to compute width of
 * \return width of tree
 */
int width(const igraph_t *tree);

/** Find the maximum delta width of a tree.
 *
 * This is the absolute value of the maximum difference in number of
 * nodes between two levels of a tree. See \cite
 * colijn2014phylogenetic.
 *
 * \param[in] tree tree to compute maximum width delta of
 * \return maximum difference between level sizes of the tree
 */
int max_delta_width(const igraph_t *tree);

/** Find the number of cherries in a tree.
 *
 * This is the number of internal nodes whose children are both
 * leaves.
 *
 * \param[in] tree tree to compute number of cherries for
 * \return number of cherries in the tree
 */
int cherries(const igraph_t *tree);

/** Find the proportion of unbalanced subtrees in a tree.
 *
 * See \cite colijn2014phylogenetic and \cite norstrom2012phylotempo.
 *
 * \param[in] tree tree to operate on
 * \return proportion of unbalanced subtrees in the tree
 */
double prop_unbalanced(const igraph_t *tree);

/** Find the average unbalance ratio of subtrees in a tree.
 *
 * This is the average of min(L, R) / max(L, R), where L (resp. R) is
 * the number of tips descending from the left (resp. right) subtree.
 * See \cite colijn2014phylogenetic and \cite norstrom2012phylotempo.
 *
 */
double avg_unbalance(const igraph_t *tree);

/** Compute Pybus' gamma statistic.
 *
 * See \cite pybus2000testing. The tree should be ultrametric, but this isn't
 * checked.
 *
 * \param[in] tree tree to compute statistic for, should be ultrametric
 * \return gamma statistic
 */
double pybus_gamma(const igraph_t *tree);

/** Compute the ratio of internal to terminal branch lengths.
 *
 * This is the ratio of the average internal branch length to the average
 * terminal branch length.
 *
 * \param[in] tree tree to compute the ratio for
 * \return ratio of internal to terminal branch lengths
 */
double internal_terminal_ratio(const igraph_t *tree);

#endif
