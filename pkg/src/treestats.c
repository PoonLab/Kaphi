#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <Judy.h>
#include <limits.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>

#include <igraph.h>

#include "tree.h"
#include "util.h"
#include "treestats.h"

#define NDEBUG

int *production(const igraph_t *tree);
int count_node_pairs(const int *production1, const int *production2, int nnode1, 
        int nnode2);
int *get_node_pairs(const igraph_t *t1, const igraph_t *t2,
        int *production1, int *production2, int npairs);
int *children(const igraph_t *tree);
double *branch_lengths(const igraph_t *tree);
double Lp_norm(const double *x1, const double *x2, const double *y1, 
        const double *y2, int n1, int n2, double p);
int _colless(const igraph_t *tree, int *n, igraph_vector_t *work, int root);
double _cophenetic(const igraph_t *tree, int *n, igraph_vector_t *work, int root, double depth, int use_branch_lengths);
int _ladder_length(const igraph_t *tree, igraph_vector_t *work, int root);
int _unbalanced(const igraph_t *tree, int *n, igraph_vector_t *work, int root);
double _unbalance(const igraph_t *tree, int *n, igraph_vector_t *work, int root);


double kernel(
    const igraph_t *t1,  // first tree as igraph object
    const igraph_t *t2,  // second tree as igraph object
    double decay_factor,  // decay constant to avoid large kernel scores when comparing tree to itself
    double rbf_variance,  // for comparing branch lengths, radial basis function variance parameter
    double sst_control,  // Moschitti's sigma (0 = subtree kernel, 1 = subset tree kernel)
    const long int *label1,  // vector of node labels converted into integers by get_labels()
    const long int *label2,
    double label_factor  // 0 completely penalizes label mismatch
    )
{

    fprintf (stdout, "entered kernel()\n");

    int i,  // handy dandy counter
        l1, l2,  //  child node indices
        r1, r2,
        n1, n2,  // internal node indices
        npairs,  // total number of pairwise node comparisons between trees
        do_label = label1 != NULL && label2 != NULL;

    int *production1 = production(t1),  // integer vector of production types
        *production2 = production(t2),
        *children1 = children(t1),  // integer indices to child nodes (twice length of num nodes)
        *children2 = children(t2);

    fprintf (stdout, "productions and children declared\n");

    double val, val1, val2, tmp, K = 0;
    double *bl1 = branch_lengths(t1),  // same length as children vector (two entries per internal node)
           *bl2 = branch_lengths(t2);

    fprintf (stdout, "branch_lengths declared\n");

    Pvoid_t delta = (Pvoid_t) NULL;  // this is a Judy array
    int *pairs, *map, cur = 0;
    PWord_t Pvalue;  // pointer into Judy array
    Word_t bytes = 0;

    fprintf (stdout, "past variable declarations\n");

    // preconditions
    assert(decay_factor > 0.0 && decay_factor <= 1.0);
    assert(rbf_variance > 0.0);
    assert(igraph_vcount(t1) < 65535);  // TODO: relax this constraint

    npairs = count_node_pairs(production1, production2, igraph_vcount(t1), igraph_vcount(t2));
    pairs = get_node_pairs(t1, t2, production1, production2, npairs);

    fprintf (stdout, "before main loop\n");

    for (cur = 0; cur < npairs; ++cur)
    {
        val = decay_factor;

        // retrieve indices for this node pair
        n1 = pairs[cur] >> 16;  // shift right by 16 bits - why?
        n2 = pairs[cur] & 65535;  // make sure we're in left side of indices

        // pre-compute branch length penalty
        tmp = pow(bl1[2*n1] - bl2[2*n2], 2) + pow(bl1[2*n1+1] - bl2[2*n2+1], 2);

        // first go down left branches and then right branches
        l1 = children1[2*n1]; r1 = children1[2*n1+1];
        l2 = children2[2*n2]; r2 = children2[2*n2+1];

        // handle paired cherries
        if (production1[l1]==0 && production1[r1]==0
            && production2[l2]==0 && production2[r2]==0) {
            // evaluate two possible comparisons of cherries
            val1 = exp(-tmp/rbf_variance);
            tmp = pow(bl1[2*n1] - bl2[2*n2+1], 2) + pow(bl1[2*n1+1] - bl2[2*n2], 2);
            val2 = exp(-tmp/rbf_variance);

            if (do_label) {
                if (label1[l1]==label2[l2] && label1[r1]==label2[r2]) {
                    // evaluating both left and right
                    val1 *= (sst_control + decay_factor) * (sst_control + decay_factor);
                } else {
                    // labels don't match
                    val1 *= (sst_control + decay_factor * label_factor) * (sst_control + decay_factor * label_factor);
                }
                if (label1[l1]==label2[r2] && label1[r1]==label2[l2]) {
                    val2 *= (sst_control + decay_factor) * (sst_control + decay_factor);
                } else {
                    // rotated labels don't match
                    val2 *= (sst_control + decay_factor * label_factor) * (sst_control + decay_factor * label_factor);
                }
            } else {
                val1 *= (sst_control + decay_factor) * (sst_control + decay_factor);
                val2 *= (sst_control + decay_factor) * (sst_control + decay_factor);
            }

            // which maximizes score?
            val *= (val1 > val2) ? val1 : val2;
        } else {
            // branch lengths
            val *= exp(-tmp/rbf_variance);

            if (production1[l1] == production2[l2]) {  // do lefts
                if (production1[l1] == 0) {  // children are leaves
                    if (do_label && label1[l1] != label2[l2]) val *= (sst_control + decay_factor * label_factor);
                    else val *= (sst_control + decay_factor);
                } else {
                    JLG(Pvalue, delta, (l1 << 16) | l2);  // get pointer Pvalue associated with index in Judy array delta
                    assert(Pvalue != NULL);  // don't visit parents before children
                    memcpy(&tmp, Pvalue, sizeof(double));
                    val *= (sst_control + tmp);
                }
            }
            if (production1[r1] == production2[r2]) {  // do rights
                if (production1[r1] == 0) {
                    if (do_label && label1[r1] != label2[r2]) val *= (sst_control + decay_factor * label_factor);
                    else val *= (sst_control + decay_factor);
                } else {
                    JLG(Pvalue, delta, (r1 << 16) | r2);
                    assert(Pvalue != NULL);
                    memcpy(&tmp, Pvalue, sizeof(double));
                    val *= (sst_control + tmp);
                }
            }
        }

        // insert into Judy array
        JLI(Pvalue, delta, pairs[cur]);
        if (Pvalue == PJERR) exit(EXIT_FAILURE);  // malloc fail occurred
        memcpy(Pvalue, &val, sizeof(double));

        K += val;
    }

    fprintf (stdout, "out of main loop\n");

    free(production1);
    free(production2);
    free(children1);
    free(children2);
    free(bl1);
    free(bl2);
    free(pairs);
    JLFA(bytes, delta);  // free the entire Judy array
    return K;
}


/*
  nLTT - normalized lineages through time measure

  Toni, T., Welch, D., Strelkowa, N., Ipsen, A., & Stumpf, M.P.H. (2009). Approximate Bayesian
  computation scheme for parameter inference and model selection in dynamical systems. Journal of
  the Royal Society Interface, 6(31), 187-202.

  We have modified this measure from the original so that it should be able to handle trees
  that are not ultrametric.
 */
double nLTT(const igraph_t *t1, const igraph_t *t2)
{
    // TODO: this needs more documentation!
    int itree, i, cur;
    double h, k, prev;
    const igraph_t *trees[2] = {t1, t2};

    // number of tips
    int n[2] = {(igraph_vcount(t1) + 1) / 2, (igraph_vcount(t2) + 1) / 2};

    // x-axis of LTT plot is node depth (distance from root in branch length units)
    double *x[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};

    // y-axis of LTT plot is number of lineages
    double *y[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};

    // array to store node depths
    double *buf = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(double));

    // to store index permutation that sorts nodes by depth
    int *node_order = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(int));

    igraph_vector_t vec;
    igraph_vector_init(&vec, igraph_vcount(t1));  // what about t2?

    for (itree = 0; itree < 2; ++itree)
    {
        depths(trees[itree], 1, buf);  // node distances from root in units of branch length
        order(buf, node_order, sizeof(double), igraph_vcount(trees[itree]),
                compare_doubles);
        // store the outdegree of the node
        igraph_degree(trees[itree], &vec, igraph_vss_all(), IGRAPH_OUT, 0);

        prev = 0; cur = 0; h = 0;
        y[itree][0] = -1.0 / (n[itree] - 1);  // anticipates increment below (root y=0)

        // iterate over every node in tree
        for (i = 0; i < igraph_vcount(trees[itree]); ++i)
        {
            // evaluate only internal nodes (outdegree > 0)
            if (VECTOR(vec)[node_order[i]] > 0) {
                if (buf[node_order[i]] == prev) {
                    // current node has same depth as previous, don't update x
                    y[itree][cur] += 1.0 / (n[itree] - 1);
                }
                else {
                    x[itree][++cur] = buf[node_order[i]];
                    y[itree][cur] = y[itree][cur-1] + 1.0 / (n[itree] - 1);
                    h = fmax(h, buf[node_order[i]]);
                    prev = buf[node_order[i]];
                }
            }
        }

        // normalize x-axis by maximum node depth/height
        for (i = 0; i < n[itree]; ++i) {  // shouldn't this be igraph_vcount, not [n]?
            x[itree][i] /= h;
        }

        // debugging
        for (i = 0; i < n[itree]-1; i++) {
            fprintf(stdout, "%d\t%d\t%d\t%f\t%f\n", itree, i, buf[node_order[i]], x[itree][i], y[itree][i]);
        }
    }

    k = Lp_norm(x[0], x[1], y[0], y[1], n[0], n[1], 1.0);

    free(x[0]);
    free(x[1]);
    free(y[0]);
    free(y[1]);
    free(buf);
    free(node_order);
    igraph_vector_destroy(&vec);
    return k;
}

double sackin(const igraph_t *t, int use_branch_lengths)
{
    int i, ntip = (igraph_vcount(t) + 1) / 2;
    double s = 0, h;
    double *buf = malloc(igraph_vcount(t) * sizeof(double));
    igraph_vector_t out_degree;

    igraph_vector_init(&out_degree, igraph_vcount(t));
    depths(t, use_branch_lengths, buf);
    igraph_degree(t, &out_degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_vcount(t); ++i) {
        if (VECTOR(out_degree)[i] == 0) {
            s += buf[i];
        }
    }

    free(buf);
    igraph_vector_destroy(&out_degree);
    return s;
}

int colless(const igraph_t *t)
{
    int c;
    int ntip = (igraph_vcount(t) + 1) / 2;
    int *n = malloc(igraph_vcount(t) * sizeof(int));
    igraph_vector_t work;
    igraph_vector_init(&work, 0);

    c = _colless(t, n, &work, root(t));
    
    igraph_vector_destroy(&work);
    free(n);
    return c;
}

double cophenetic(const igraph_t *tree, int use_branch_lengths)
{
    double c, h;
    igraph_vector_t work;
    igraph_vector_init(&work, 0);
    int *n = malloc(igraph_vcount(tree) * sizeof(int));
    int ntip = (igraph_vcount(tree) + 1) / 2;

    c = _cophenetic(tree, n, &work, root(tree), 0, use_branch_lengths);

    free(n);
    igraph_vector_destroy(&work);
    return c;
}

int ladder_length(const igraph_t *tree)
{
    int l;
    igraph_vector_t work;
    igraph_vector_init(&work, 0);

    l = _ladder_length(tree, &work, root(tree));

    igraph_vector_destroy(&work);
    return l;
}

int il_nodes(const igraph_t *tree)
{
    int i, lc, rc;
    int il = 0;
    igraph_vector_t degree, child;

    igraph_vector_init(&degree, igraph_vcount(tree));
    igraph_vector_init(&child, 2);
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        if (VECTOR(degree)[i] > 0) {
            igraph_neighbors(tree, &child, i, IGRAPH_OUT);
            lc = VECTOR(child)[0]; rc = VECTOR(child)[1];
            il += VECTOR(degree)[lc] == 0 ^ VECTOR(degree)[rc] == 0;
        }
    }

    igraph_vector_destroy(&child);
    igraph_vector_destroy(&degree);
    return il;
}

int width(const igraph_t *tree)
{
    int i = 0, w = 0, wprev = 1;
    igraph_vs_t vs = igraph_vss_1(root(tree));

    igraph_vector_t nbhd;
    igraph_vector_init(&nbhd, 0);
    igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT);

    while (VECTOR(nbhd)[0] < igraph_vcount(tree))
    {
        ++i;
        igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT);
        if (VECTOR(nbhd)[0] - wprev > w) {
            w = VECTOR(nbhd)[0] - wprev;
        }
        wprev = VECTOR(nbhd)[0];
    }
    igraph_vector_destroy(&nbhd);
    return w;
} 

int max_delta_width(const igraph_t *tree)
{
    int i = 0, w = 0, wprev = 1, nprev = 1, maxdw = 0;
    igraph_vs_t vs = igraph_vss_1(root(tree));

    igraph_vector_t nbhd;
    igraph_vector_init(&nbhd, 0);
    igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT);

    while (VECTOR(nbhd)[0] < igraph_vcount(tree))
    {
        igraph_neighborhood_size(tree, &nbhd, vs, ++i, IGRAPH_OUT);
        w = VECTOR(nbhd)[0] - nprev;
        if (fabs(w - wprev) > maxdw) {
            maxdw = w > wprev ? w - wprev : wprev - w;
        }
        wprev = w;
        nprev = VECTOR(nbhd)[0];
    }
    igraph_vector_destroy(&nbhd);
    return maxdw;
}

int cherries(const igraph_t *tree)
{
    int i, lc, rc;
    int c = 0;
    igraph_vector_t degree, child;

    igraph_vector_init(&degree, igraph_vcount(tree));
    igraph_vector_init(&child, 2);
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        if (VECTOR(degree)[i] > 0) {
            igraph_neighbors(tree, &child, i, IGRAPH_OUT);
            lc = VECTOR(child)[0]; rc = VECTOR(child)[1];
            c += VECTOR(degree)[lc] == 0 && VECTOR(degree)[rc] == 0;
        }
    }

    igraph_vector_destroy(&child);
    igraph_vector_destroy(&degree);
    return c;
}

double prop_unbalanced(const igraph_t *tree)
{
    int p;
    int nnode = (igraph_vcount(tree) - 1) / 2;
    int *n = malloc(igraph_vcount(tree) * sizeof(int));
    igraph_vector_t work;
    igraph_vector_init(&work, 0);

    p = _unbalanced(tree, n, &work, root(tree));
    
    igraph_vector_destroy(&work);
    free(n);
    return (double) p / (double) nnode;
}

double avg_unbalance(const igraph_t *tree)
{
    double u;
    double nnode = (igraph_vcount(tree) - 1) / 2;
    int *n = malloc(igraph_vcount(tree) * sizeof(int));
    igraph_vector_t work;
    igraph_vector_init(&work, 0);

    u = _unbalance(tree, n, &work, root(tree));
    
    igraph_vector_destroy(&work);
    free(n);
    return u / nnode;
}

double pybus_gamma(const igraph_t *tree)
{
    int i, j, k;
    double *node_depths = malloc(igraph_vcount(tree) * sizeof(double));
    double *T = calloc((NTIP(tree) - 1), sizeof(double));
    double prev_depth = 0, gamma = 0;
    int *node_order = malloc(igraph_vcount(tree) * sizeof(int));
    igraph_vector_t degree;

    igraph_vector_init(&degree, 0);
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    depths(tree, 1, node_depths);
    order(node_depths, node_order, sizeof(double), igraph_vcount(tree), compare_doubles);

    // T[k] = sum_(j=2)^(k-2) j*g_j
    j = 2;
    for (i = 1; i < igraph_vcount(tree); ++i)
    {
        if (VECTOR(degree)[node_order[i]] > 0) {
            for (k = j - 2; k < (NTIP(tree) - 1); ++k)
            {
                T[k] += j * (node_depths[node_order[i]] - prev_depth);
            }
            ++j;
            prev_depth = node_depths[node_order[i]];
        }
    }
    T[NTIP(tree) - 2] += NTIP(tree) * (height(tree) - prev_depth);

    for (i = 2; i < NTIP(tree); ++i) {
        gamma += T[i - 2];
    }
    gamma *= 1.0 / (NTIP(tree) - 2.0);
    gamma -= T[NTIP(tree) - 2] / 2.0;
    gamma /= T[NTIP(tree) - 2] * sqrt(1.0 / 12.0 / (NTIP(tree) - 2));

    igraph_vector_destroy(&degree);
    free(node_depths);
    free(T);
    free(node_order);
    return gamma;
}

double internal_terminal_ratio(const igraph_t *tree)
{
    int i, from, to;
    int ntip = (igraph_vcount(tree) + 1) / 2;
    double internal = 0, terminal = 0;
    igraph_vector_t degree;

    igraph_vector_init(&degree, 0);
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_ecount(tree); ++i) {
        igraph_edge(tree, i, &from, &to);
        if (VECTOR(degree)[to] == 0) {
            terminal += EAN(tree, "length", i);
        }
        else {
            internal += EAN(tree, "length", i);
        }
    }
    internal /= ntip - 2;
    terminal /= ntip;
    igraph_vector_destroy(&degree);
    return internal / terminal;
}

/* Private. */

/* recursively compute the ladder length */
int _ladder_length(const igraph_t *tree, igraph_vector_t *work, int root)
{
    int lc, rc; 
    int llad, rlad;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0) {
        lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
        llad = _ladder_length(tree, work, lc);
        rlad = _ladder_length(tree, work, rc);
        return llad > rlad ? llad + 1 : rlad + 1;
    }
    else {
        return 1;
    }
}

/* recursively compute unbalance */
double _unbalance(const igraph_t *tree, int *n, igraph_vector_t *work, int root)
{
    int lc, rc; 
    double lunb, runb;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0) {
        lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
        lunb = _unbalance(tree, n, work, lc);
        runb = _unbalance(tree, n, work, rc);
        n[root] = n[lc] + n[rc];
        return lunb + runb + fmin(n[lc], n[rc]) / fmax(n[lc], n[rc]);
    }
    else {
        n[root] = 1;
        return 0;
    }
}

/* recursively compute number of unbalanced subtrees */
int _unbalanced(const igraph_t *tree, int *n, igraph_vector_t *work, int root)
{
    int lc, rc; 
    int lunb, runb;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0) {
        lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
        lunb = _unbalanced(tree, n, work, lc);
        runb = _unbalanced(tree, n, work, rc);
        n[root] = n[lc] + n[rc];
        return lunb + runb + (n[lc] != n[rc]);
    }
    else {
        n[root] = 1;
        return 0;
    }
}

/* recursively compute Colless' index */
int _colless(const igraph_t *tree, int *n, igraph_vector_t *work, int root)
{
    int lc, rc; 
    int lcol, rcol;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0) {
        lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
        lcol = _colless(tree, n, work, lc);
        rcol = _colless(tree, n, work, rc);
        n[root] = n[lc] + n[rc];
        return lcol + rcol + fabs(n[lc] - n[rc]);
    }
    else {
        n[root] = 1;
        return 0;
    }
}

/* recursively compute the total cophenetic index */
double _cophenetic(const igraph_t *tree, int *n, igraph_vector_t *work, int root, double depth, int use_branch_lengths)
{
    int lc, rc, from; 
    double lphi, rphi, lbl = 1, rbl = 1;

    igraph_incident(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0) {
        igraph_edge(tree, VECTOR(*work)[0], &from, &lc);
        igraph_edge(tree, VECTOR(*work)[1], &from, &rc);
        
        if (use_branch_lengths) {
            lbl = EAN(tree, "length", VECTOR(*work)[0]);
            rbl = EAN(tree, "length", VECTOR(*work)[1]);
        }

        lphi = _cophenetic(tree, n, work, lc, depth + lbl, use_branch_lengths);
        rphi = _cophenetic(tree, n, work, rc, depth + rbl, use_branch_lengths);
        n[root] = n[lc] + n[rc];
        return lphi + rphi + n[lc] * n[rc] * depth;
    }
    else {
        n[root] = 1;
        return 0;
    }
}

double Lp_norm(const double *x1, const double *x2, const double *y1, 
        const double *y2, int n1, int n2, double p)
{
    double a, b, fa1, fa2, fb1, fb2, m1, m2, intersect, fintersect, area, norm = 0;
    int i1 = 0, i2 = 0;

    if (x1[0] < x2[0])
    {
        b = x1[0];
        fb1 = y1[0];
        fb2 = ((y2[1] - y2[0]) / (x2[1] - x2[0])) * (b - x2[0]);
    }
    else if (x1[0] == x2[0])
    {
        b = x1[0];
        fb1 = y1[0];
        fb2 = y2[0];
    }
    else
    {
        b = x2[0];
        fb2 = y2[0];
        fb1 = ((y1[1] - y1[0]) / (x1[1] - x1[0])) * (b - x1[0]);
    }

    while (i1 < n1 && i2 < n2)
    {
        a = b;
        b = fmin(x1[i1], x2[i2]);
        fa1 = fb1;
        fa2 = fb2;

        if (x1[i1] < x2[i2])
        {
            fb1 = y1[i1++];
            fb2 = fa2 + (y2[i2] - fa2) / (x2[i2] - a) * (b-a);
        }
        else if (x1[i1] == x2[i2])
        {
            fb1 = y1[i1++];
            fb2 = y2[i2++];
        }
        else
        {
            fb2 = y2[i2++];
            fb1 = fa1 + (y1[i1] - fa1) / (x1[i1] - a) * (b-a);
        }

        m1 = (fb1 - fa1) / (b - a);
        m2 = (fb2 - fa2) / (b - a);

        if (m1 != m2)
        {
            intersect = a + (fa1 - fa2) / (m2 - m1);
        }
        else
        {
            intersect = a;
        }

        if (intersect > a && intersect < b)
        {
            fintersect = fa1 + m1 * intersect;
            area  = fabs(((fa1 + fintersect) / 2 - (fa2 + fintersect) / 2) * (intersect-a));
            area += fabs(((fintersect + fb1) / 2 - (fintersect + fb2) / 2) * (b-intersect));
        }
        else
        {
            area = fabs(((fa1 + fb1) / 2 - (fa2 + fb2) / 2) * (b-a));
        }
        norm += area;
    }
    return pow(norm, p);
}


/*
    production()
    Get production rules for each node in a binary tree.
    Returns rules as an integer vector.
     0 = terminal node
     1 = internal node, both children are internal nodes
     2 = internal node, one child is a terminal node
     3 = internal node, both children are terminal nodes
*/
int *production(const igraph_t *tree)
{
    int i, nnode = igraph_vcount(tree);
    int *p = malloc(nnode * sizeof(int));
    igraph_vector_t vec, nbr;
    igraph_vector_init(&nbr, 2);
    igraph_vector_init(&vec, 2);

    igraph_degree(tree, &vec, igraph_vss_all(), IGRAPH_OUT, 0);
    for (i = 0; i < nnode; ++i)
    {
        if ((int) VECTOR(vec)[i] == 0)
        {
            // node is terminal with no out edges (a tip)
            p[i] = 0;
        }
        else
        {
            // node is internal, production index is at least 1
            igraph_neighbors(tree, &nbr, i, IGRAPH_OUT);

            // adjust index by number of child nodes that are tips
            // if neither child is a tip, index is 1
            // if one child is a tip, index is 2
            // if both children are tips, index is 3
            p[i] = ((int) VECTOR(vec)[(int) VECTOR(nbr)[0]] == 0) +
                   ((int) VECTOR(vec)[(int) VECTOR(nbr)[1]] == 0) + 1;
        }
    }
    igraph_vector_destroy(&nbr);
    igraph_vector_destroy(&vec);
    return p;
}

/* find how many pairs of nodes we need to evaluate */
int count_node_pairs(const int *production1, const int *production2, int nnode1, 
        int nnode2)
{
    int i, npairs = 0, table1[4] = {0}, table2[4] = {0};

    for (i = 0; i < nnode1; ++i) {
        table1[production1[i]] += 1;
    }
    for (i = 0; i < nnode2; ++i) {
        table2[production2[i]] += 1;
    }

    for (i = 1; i < 4; ++i) {
        npairs += table1[i] * table2[i];
    }
    return npairs;
}

/* get all pairs of nodes with the same production */
int *get_node_pairs(const igraph_t *t1, const igraph_t *t2,
        int *production1, int *production2, int npairs)
{
    int coord, n1, n2, i1 = 1, i2 = 1, start = 0;
    int nnode1 = igraph_vcount(t1), nnode2 = igraph_vcount(t2);
    int *order1, *order2;
    int *pairs, i = 0;

    order1 = malloc(nnode1 * sizeof(int));
    order2 = malloc(nnode2 * sizeof(int));
    order(production1, order1, sizeof(int), nnode1, compare_ints);
    order(production2, order2, sizeof(int), nnode2, compare_ints);

    pairs = malloc(npairs * sizeof(int));

    n1 = order1[0]; n2 = order2[0];

    // find pairs of nodes with equal productions
    while (i1 < nnode1 || i2 < nnode2) 
    {
        while (production1[n1] != production2[n2])
        {
            if (production1[n1] > production2[n2])
            {
                if (i2 < nnode2)
                    n2 = order2[i2++];
                else
                    break;
            }
            else if (production1[n1] < production2[n2])
            {
                if (i1 < nnode1)
                    n1 = order1[i1++];
                else
                    break;
            }
        } 
        start = i2 - 1;

        if (production1[n1] != production2[n2])
            break;

        while (production1[n1] == production2[n2])
        {
            while (production1[n1] == production2[n2])
            {
                coord = (n1 << 16) | n2;

                // internal nodes
                if (production1[n1] > 0) 
                    pairs[i++] = coord;

                if (i2 < nnode2) 
                    n2 = order2[i2++];
                else 
                    break;
            }
            if (i1 < nnode1)
            {
                n1 = order1[i1++];
                i2 = start;
                n2 = order2[i2++];
            }
            else
            {
                break;
            }
        }
    }
    free(order1);
    free(order2);

    qsort(pairs, npairs, sizeof(int), compare_ints);
    return pairs;
}

/* get indices of children from each node */
int *children(const igraph_t *tree)
{
    igraph_adjlist_t al;
    igraph_vector_int_t *nbr;
    int i;
    int *children = malloc(igraph_vcount(tree) * 2 * sizeof(int));

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        nbr = igraph_adjlist_get(&al, i);
        if (igraph_vector_int_size(nbr) > 0)
        {
            children[2*i] = VECTOR(*nbr)[0];
            children[2*i+1] = VECTOR(*nbr)[1];
        }
    }

    igraph_adjlist_destroy(&al);
    return children;
}

/* get branch lengths leading out of each node */
double *branch_lengths(const igraph_t *tree)
{
    fprintf (stdout, "entered branch_lengths()\n");

    igraph_inclist_t il;  // incidence list
    int i;
    igraph_vector_int_t *edge;
    double *branch_lengths = malloc(2 * igraph_vcount(tree) * sizeof(double));

    fprintf (stdout, "past variable declarations\n");

    igraph_inclist_init(tree, &il, IGRAPH_OUT);

    fprintf (stdout, "past inclist_init, about to enter main loop\n");
    fprintf (stdout, "vcount: %d\n", igraph_vcount(tree));

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        fprintf (stdout, "%d\n", i);

        edge = igraph_inclist_get(&il, i);
        if (igraph_vector_int_size(edge) > 0)
        {
            branch_lengths[2*i] = EAN(tree, "length", VECTOR(*edge)[0]);
            fprintf (stdout, " left %f\n", branch_lengths[2*i]);

            branch_lengths[2*i+1] = EAN(tree, "length", VECTOR(*edge)[1]);
            fprintf (stdout, )
            fprintf (stdout, " right %f\n", branch_lengths[2*i+1]);
        }
    }

    fprintf (stdout, "exited main loop\n");

    igraph_inclist_destroy(&il);

    fprintf (stdout, "past inclist_destroy(), returning from branch_lengths()\n");

    return branch_lengths;
}




/* from en.wikipedia.org/wiki/https://en.wikipedia.org/wiki/Heapsort */
void my_igraph_strvector_swap(igraph_strvector_t *a, long int i1, long int i2) {
    char *tmp1;
    char *tmp2;
    char *tmp3;
    char *tmp4;

    igraph_strvector_get(a, i1, &tmp1);
    igraph_strvector_get(a, i2, &tmp2);
    tmp3 = (char*)malloc(sizeof(char*) * (strlen(tmp1) + 1));
    strcpy(tmp3, tmp1);
    tmp4 = (char*)malloc(sizeof(char*) * (strlen(tmp2) + 1));
    strcpy(tmp4, tmp2);
    igraph_strvector_set(a, i1, tmp4);
    igraph_strvector_set(a, i2, tmp3);
    free(tmp3);
    free(tmp4);
}

void sift_down(igraph_strvector_t *a, long int start, long int end) {
    long int root = start;
    long int child, swap;
    char *tmp1;
    char *tmp2;

    while (2 * root + 1 <= end) {
        child = 2 *root + 1;
        swap = root;

        igraph_strvector_get(a, swap, &tmp1);
        igraph_strvector_get(a, child, &tmp2);
        if (strcmp(tmp1, tmp2) < 0) {
            swap = child;
        }

        if (child + 1 <= end) {
            igraph_strvector_get(a, swap, &tmp1);
            igraph_strvector_get(a, child + 1, &tmp2);
            if (strcmp(tmp1, tmp2) < 0) {
                swap = child + 1;
            }
        }

        if (swap == root)
            return;
        else {
            my_igraph_strvector_swap(a, swap, root);

            root = swap;
        }
    }
}

void heapify(igraph_strvector_t *a, long int count) {
    long int start = (count - 2) / 2;

    while (start >= 0L) {
        sift_down(a, start, count - 1);
        start = start - 1;
    }
}

void my_igraph_strvector_sort(igraph_strvector_t *a) {
    long int count = igraph_strvector_size(a);
    long int end;

    heapify(a, count);
    end = count - 1;

    while (end > 0L) {
        my_igraph_strvector_swap(a, 0, end);

        end = end - 1;
        sift_down(a, 0L, end);
    }
}

long int my_igraph_strvector_search(const igraph_strvector_t *a, char *val) {
    long int start = 0;
    long int end = igraph_strvector_size(a);
    long int mid = (start + end) / 2;
    int cmp;
    char *tmp;

    while (start < end) {
        igraph_strvector_get(a, mid, &tmp);

        cmp = strcmp(val, tmp);
        if (cmp < 0) {
            end = mid;
        }
        else {
            if (cmp > 0) {
                start = mid + 1;
            }
            else {
                // cmp == 0, both strings are equal
                return mid;
            }
        }

        mid = (start + end) / 2;
    }

    return -1;
}

/* convert tree labels to index reference */
void get_labels(const igraph_t *tree1, const igraph_t *tree2, long int *label1, long int *label2) {
    long int i;
    char * current_string;

    // initialize containers for node labels
    igraph_strvector_t string_label1;
    igraph_strvector_t string_label2;
    igraph_strvector_t string_label_all;

    igraph_strvector_init(&string_label1, 0);
    igraph_strvector_init(&string_label2, 0);

    // igraph - query string vertex attribute for all vertices
    VASV(tree1, "id", &string_label1);
    VASV(tree2, "id", &string_label2);

    // concatenate string vectors
    igraph_strvector_init(&string_label_all, 0);
    igraph_strvector_copy(&string_label_all, &string_label1);
    igraph_strvector_append(&string_label_all, &string_label2);

    // heap sort of all node labels
    my_igraph_strvector_sort(&string_label_all);

    // convert node labels into a reduced set of integer values based on heap sort
    for (i = 0; i < igraph_strvector_size(&string_label1); ++i) {
        igraph_strvector_get(&string_label1, i, &current_string);
        label1[i] = my_igraph_strvector_search(&string_label_all, current_string);
        //printf("%ld: '%s' (%ld)\n", i, current_string, label1[i]);
    }

    for (i = 0; i < igraph_strvector_size(&string_label2); ++i) {
        igraph_strvector_get(&string_label2, i, &current_string);
        label2[i] = my_igraph_strvector_search(&string_label_all, current_string);
        //printf("%ld: '%s' (%ld)\n", i, current_string, label2[i]);
    }

    //printf("\n");

    // free memory
    igraph_strvector_destroy(&string_label_all);
    igraph_strvector_destroy(&string_label1);
    igraph_strvector_destroy(&string_label2);
}


