#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <Judy.h>
#include <limits.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>

#include "../igraph/include/igraph.h"

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

double kernel(const igraph_t *t1, const igraph_t *t2, double decay_factor, 
        double rbf_variance, double sst_control)
{
    int i, c1, c2, n1, n2, coord, npairs;
    int *production1, *production2, *children1, *children2;
    double val, tmp, K = 0;
    double *bl1, *bl2;
    Pvoid_t delta = (Pvoid_t) NULL;
    int *pairs, *map, cur = 0;
    PWord_t Pvalue;
    Word_t bytes = 0;

    // preconditions
    assert(decay_factor > 0.0 && decay_factor <= 1.0);
    assert(rbf_variance > 0.0);
    assert(igraph_vcount(t1) < 65535);

    production1 = production(t1);
    production2 = production(t2);
    children1 = children(t1);
    children2 = children(t2);
    bl1 = branch_lengths(t1);
    bl2 = branch_lengths(t2);

    npairs = count_node_pairs(production1, production2, igraph_vcount(t1), igraph_vcount(t2));
    pairs = get_node_pairs(t1, t2, production1, production2, npairs);
    
    for (cur = 0; cur < npairs; ++cur)
    {
        val = decay_factor;
        n1 = pairs[cur] >> 16;
        n2 = pairs[cur] & 65535;

        // branch lengths
        tmp = pow(bl1[2*n1] - bl2[2*n2], 2) + pow(bl1[2*n1+1] - bl2[2*n2+1], 2);
        val *= exp(-tmp/rbf_variance);

        for (i = 0; i < 2; ++i)  // assume tree is binary
        {
            c1 = children1[2*n1+i];
            c2 = children2[2*n2+i];

            if (production1[c1] == production2[c2])
            {
                // children are leaves
                if (production1[c1] == 0)
                {
                    val *= (sst_control + decay_factor);
                }

                // children are not leaves
                else
                {
                    JLG(Pvalue, delta, (c1 << 16) | c2);
                    /* don't visit parents before children */
                    assert(Pvalue != NULL);
                    memcpy(&tmp, Pvalue, sizeof(double));
                    val *= (sst_control + tmp);
                }
            }
        }

        JLI(Pvalue, delta, pairs[cur]);
        if (Pvalue == PJERR) exit(EXIT_FAILURE); 
        memcpy(Pvalue, &val, sizeof(double));

        K += val;
    }

    free(production1);
    free(production2);
    free(children1);
    free(children2);
    free(bl1);
    free(bl2);
    free(pairs);
    JLFA(bytes, delta);
    return K;
}

double nLTT(const igraph_t *t1, const igraph_t *t2)
{
    int itree, i, cur;
    double h, k, prev;
    const igraph_t *trees[2] = {t1, t2};

    int n[2] = {(igraph_vcount(t1) + 1) / 2, (igraph_vcount(t2) + 1) / 2};
    double *x[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};
    double *y[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};
    double *buf = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(double));
    int *node_order = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(int));

    igraph_vector_t vec;
    igraph_vector_init(&vec, igraph_vcount(t1));

    for (itree = 0; itree < 2; ++itree)
    {
        depths(trees[itree], 1, buf);
        order(buf, node_order, sizeof(double), igraph_vcount(trees[itree]),
                compare_doubles);
        igraph_degree(trees[itree], &vec, igraph_vss_all(), IGRAPH_OUT, 0);

        prev = 0; cur = 0; h = 0;
        y[itree][0] = -1.0 / (n[itree] - 2);
        for (i = 0; i < igraph_vcount(trees[itree]); ++i)
        {
            if (VECTOR(vec)[node_order[i]] > 0) {
                if (buf[node_order[i]] == prev) {
                    y[itree][cur] += 1.0 / (n[itree] - 2);
                }
                else {
                    x[itree][++cur] = buf[node_order[i]];
                    y[itree][cur] = y[itree][cur-1] + 1.0 / (n[itree] - 2);
                    h = fmax(h, buf[node_order[i]]);
                    prev = buf[node_order[i]];
                }
            }
        }

        for (i = 0; i < n[itree]; ++i) {
            x[itree][i] /= h;
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
    igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT, 0);

    while (VECTOR(nbhd)[0] < igraph_vcount(tree))
    {
        ++i;
        igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT, 0);
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
    igraph_neighborhood_size(tree, &nbhd, vs, i, IGRAPH_OUT, 0);

    while (VECTOR(nbhd)[0] < igraph_vcount(tree))
    {
        igraph_neighborhood_size(tree, &nbhd, vs, ++i, IGRAPH_OUT, 0);
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

/* get production rules for each node */
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
            p[i] = 0;
        }
        else
        {
            igraph_neighbors(tree, &nbr, i, IGRAPH_OUT);
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
    igraph_inclist_t il;
    int i;
    igraph_vector_int_t *edge;
    double *branch_lengths = malloc(2 * igraph_vcount(tree) * sizeof(double));

    igraph_inclist_init(tree, &il, IGRAPH_OUT);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        edge = igraph_inclist_get(&il, i);
        if (igraph_vector_int_size(edge) > 0)
        {
            branch_lengths[2*i] = EAN(tree, "length", VECTOR(*edge)[0]);
            branch_lengths[2*i+1] = EAN(tree, "length", VECTOR(*edge)[1]);
        }
    }

    igraph_inclist_destroy(&il);
    return branch_lengths;
}

