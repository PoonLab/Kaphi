#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "../igraph/include/igraph.h"

#include "newick_parser.h"
#include "tree.h"
#include "util.h"

#define NDEBUG

void yyrestart(FILE *f);

typedef struct {
    igraph_strvector_t *strvattrs;
    igraph_vector_t *numvattrs;
    igraph_vector_bool_t *boolvattrs;

    igraph_strvector_t *streattrs;
    igraph_vector_t *numeattrs;
    igraph_vector_bool_t *booleattrs;

    int nstrvattr;
    int nnumvattr;
    int nboolvattr;

    int nstreattr;
    int nnumeattr;
    int nbooleattr;

    igraph_strvector_t strvattr_names;
    igraph_strvector_t numvattr_names;
    igraph_strvector_t boolvattr_names;

    igraph_strvector_t streattr_names;
    igraph_strvector_t numeattr_names;
    igraph_strvector_t booleattr_names;
} tree_attrs;

/* recursive functions */
int _ladderize(igraph_t *tree, igraph_vector_t *work, int root, int *perm);
double _height(const igraph_t *tree, igraph_vector_t *work, int root);
int _write_tree_newick(const igraph_t *tree, char *out, int root, igraph_vector_t *vec);
void _cut_at_time(igraph_t *tree, double t, int root, double troot, 
        int extant_only, igraph_vector_t *work, igraph_vector_t *to_delete);
int _collapse_singles(igraph_t *tree, int root, igraph_vector_t *vdel,
        igraph_vector_t *eadd, igraph_vector_t *branch_length,
        igraph_vector_t *work, double *bl);
void _depths(const igraph_t *tree, double *depths, igraph_vector_t *work, 
        int root, double parent_depth, int use_branch_lengths);

/* other helpers */
tree_attrs *_get_tree_attrs(const igraph_t *tree);
void _permute_tree_attrs(igraph_t *tree, tree_attrs *a, const int *perm);
void _tree_attrs_destroy(tree_attrs *a);
void _get_node_ids(const igraph_t *g, igraph_strvector_t *ids);
void _make_maps(const igraph_t *tree, const igraph_t *net, int *tip_map, int *node_map);

igraph_t *parse_newick(FILE *f)
{
    igraph_t *tree;
    int i;
    extern int yynode;
    igraph_vector_t edge, branch_length, size;
    igraph_strvector_t label;

    igraph_vector_init(&edge, 0);
    igraph_vector_init(&size, 0);
    igraph_vector_init(&branch_length, 0);
    igraph_strvector_init(&label, 0);

    yynode = 0;
    yyrestart(f);
    yyparse(&edge, &size, &branch_length, &label);

    tree = malloc(sizeof(igraph_t));
    igraph_empty(tree, igraph_vector_size(&size), 1);
    igraph_add_edges(tree, &edge, 0);

    for (i = 0; i < igraph_vector_size(&size); ++i)
    {
        igraph_incident(tree, &edge, i, IGRAPH_IN);
        if (igraph_vector_size(&edge) > 0) {
            SETEAN(tree, "length", (int) VECTOR(edge)[0], VECTOR(branch_length)[i]);
        }
        SETVAS(tree, "id", i, STR(label,i));
    }

    igraph_vector_destroy(&edge);
    igraph_vector_destroy(&size);
    igraph_vector_destroy(&branch_length);
    igraph_strvector_destroy(&label);
    return tree;
}

int root(const igraph_t *tree)
{
    long r = -1;
    igraph_vector_t work;

    igraph_vector_init(&work, 1);
    igraph_degree(tree, &work, igraph_vss_all(), IGRAPH_IN, 0);
    igraph_vector_search(&work, 0, 0, &r);
    igraph_vector_destroy(&work);
    return (int) r;
}

void ladderize(igraph_t *tree)
{
    igraph_vector_t work, permvec;
    igraph_t *new_tree = malloc(sizeof(igraph_t));
    int i;
    int *order = malloc(igraph_vcount(tree) * sizeof(int));

    igraph_vector_init(&work, 2);
    igraph_vector_init(&permvec, igraph_vcount(tree));

    _ladderize(tree, &work, root(tree), order);

    for (i = 0; i < igraph_vcount(tree); ++i) {
        igraph_vector_set(&permvec, order[i], (double) i);
    }

    igraph_permute_vertices(tree, new_tree, &permvec);
    igraph_destroy(tree);
    memcpy(tree, new_tree, sizeof(igraph_t));

    free(order);
    free(new_tree);
    igraph_vector_destroy(&work);
    igraph_vector_destroy(&permvec);
}

tree_attrs *_get_tree_attrs(const igraph_t *tree)
{
    int i, j, nv = igraph_vcount(tree), ne = igraph_ecount(tree), from, to;
    int *heads = malloc(ne * sizeof(int));
    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_t gtypes, vtypes, etypes;
    tree_attrs *a = malloc(sizeof(tree_attrs));

    igraph_strvector_t streattr;
    igraph_vector_t numeattr;
    igraph_vector_bool_t booleattr;

    for (i = 0; i < ne; ++i) {
        igraph_edge(tree, i, &from, &to);
        heads[i] = to;
    }

    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, nv);
    igraph_strvector_init(&enames, nv);

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, nv);
    igraph_vector_init(&etypes, nv);

    igraph_strvector_init(&streattr, 0);
    igraph_vector_init(&numeattr, 0);
    igraph_vector_bool_init(&booleattr, 0);

    igraph_cattribute_list(tree, &gnames, &gtypes, &vnames, &vtypes, &enames,
                           &etypes);

    igraph_strvector_init(&a->strvattr_names, 0);
    igraph_strvector_init(&a->numvattr_names, 0);
    igraph_strvector_init(&a->boolvattr_names, 0);

    igraph_strvector_init(&a->streattr_names, 0);
    igraph_strvector_init(&a->numeattr_names, 0);
    igraph_strvector_init(&a->booleattr_names, 0);

    for (i = 0; i < igraph_vector_size(&vtypes); ++i) {
        switch ((igraph_attribute_type_t) VECTOR(vtypes)[i]) {
            case IGRAPH_ATTRIBUTE_STRING:
                igraph_strvector_add(&a->strvattr_names, STR(vnames, i));
                break;
            case IGRAPH_ATTRIBUTE_NUMERIC:
                igraph_strvector_add(&a->numvattr_names, STR(vnames, i));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                igraph_strvector_add(&a->boolvattr_names, STR(vnames, i));
                break;
            default:
                break;
        }
    }

    for (i = 0; i < igraph_vector_size(&etypes); ++i) {
        switch ((igraph_attribute_type_t) VECTOR(etypes)[i]) {
            case IGRAPH_ATTRIBUTE_STRING:
                igraph_strvector_add(&a->streattr_names, STR(enames, i));
                break;
            case IGRAPH_ATTRIBUTE_NUMERIC:
                igraph_strvector_add(&a->numeattr_names, STR(enames, i));
                break;
            case IGRAPH_ATTRIBUTE_BOOLEAN:
                igraph_strvector_add(&a->booleattr_names, STR(enames, i));
                break;
            default:
                break;
        }
    }

    a->nstrvattr = igraph_strvector_size(&a->strvattr_names);
    a->nnumvattr = igraph_strvector_size(&a->numvattr_names);
    a->nboolvattr = igraph_strvector_size(&a->boolvattr_names);

    a->nstreattr = igraph_strvector_size(&a->streattr_names);
    a->nnumeattr = igraph_strvector_size(&a->numeattr_names);
    a->nbooleattr = igraph_strvector_size(&a->booleattr_names);

    a->strvattrs = malloc(a->nstrvattr * sizeof(igraph_strvector_t));
    a->numvattrs = malloc(a->nnumvattr * sizeof(igraph_vector_t));
    a->boolvattrs = malloc(a->nboolvattr * sizeof(igraph_vector_bool_t));

    a->streattrs = malloc(a->nstreattr * sizeof(igraph_strvector_t));
    a->numeattrs = malloc(a->nnumeattr * sizeof(igraph_vector_t));
    a->booleattrs = malloc(a->nbooleattr * sizeof(igraph_vector_bool_t));

    for (i = 0; i < a->nstrvattr; ++i) {
        igraph_strvector_init(&a->strvattrs[i], nv);
        VASV(tree, STR(a->strvattr_names, i), &a->strvattrs[i]);
    }
    for (i = 0; i < a->nnumvattr; ++i) {
        igraph_vector_init(&a->numvattrs[i], nv);
        VANV(tree, STR(a->strvattr_names, i), &a->numvattrs[i]);
    }
    for (i = 0; i < a->nboolvattr; ++i) {
        igraph_vector_bool_init(&a->boolvattrs[i], nv);
        VABV(tree, STR(a->boolvattr_names, i), &a->boolvattrs[i]);
    }

    for (i = 0; i < a->nstreattr; ++i) {
        igraph_strvector_init(&a->streattrs[i], nv);
        EASV(tree, STR(a->streattr_names, i), &streattr);
        for (j = 0; j < ne; ++j) {
            igraph_strvector_set(&a->streattrs[i], heads[j], STR(streattr, j));
        }
    }
    for (i = 0; i < a->nnumeattr; ++i) {
        igraph_vector_init(&a->numeattrs[i], nv);
        EANV(tree, STR(a->numeattr_names, i), &numeattr);
        for (j = 0; j < ne; ++j) {
            VECTOR(a->numeattrs[i])[heads[j]] = VECTOR(numeattr)[j];
        }
    }
    for (i = 0; i < a->nbooleattr; ++i) {
        igraph_vector_bool_init(&a->booleattrs[i], nv);
        EABV(tree, STR(a->booleattr_names, i), &booleattr);
        for (j = 0; j < ne; ++j) {
            VECTOR(a->booleattrs[i])[heads[j]] = VECTOR(booleattr)[j];
        }
    }

    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);

    igraph_vector_destroy(&gtypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&etypes);

    igraph_strvector_destroy(&streattr);
    igraph_vector_destroy(&numeattr);
    igraph_vector_bool_destroy(&booleattr);

    free(heads);

    return a;
}

void _permute_tree_attrs(igraph_t *tree, tree_attrs *a, const int *perm)
{
    int i, j, nv = igraph_vcount(tree);
    igraph_vector_int_t *inc;
    igraph_inclist_t inclist;
    int *edge = malloc(nv * sizeof(int));

    igraph_inclist_init(tree, &inclist, IGRAPH_IN);
    for (i = 0; i < nv; ++i) {
        inc = igraph_inclist_get(&inclist, i);
        if (igraph_vector_int_size(inc) > 0) {
            edge[i] = VECTOR(*inc)[0];
        }
        else {
            edge[i] = -1;
        }
    }
    igraph_inclist_destroy(&inclist);

    for (i = 0; i < a->nstrvattr; ++i) {
        permute(&a->strvattrs[i], BUFSIZ, nv, perm, get_igraph_strvector_t, 
                set_igraph_strvector_t);
        SETVASV(tree, STR(a->strvattr_names, i), &a->strvattrs[i]);
    }

    for (i = 0; i < a->nnumvattr; ++i) {
        permute(&a->numvattrs[i], sizeof(igraph_real_t), nv, perm, 
                get_igraph_vector_t, set_igraph_vector_t);
        SETVANV(tree, STR(a->numvattr_names, i), &a->numvattrs[i]);
    }

    for (i = 0; i < a->nboolvattr; ++i) {
        permute(&a->boolvattrs[i], sizeof(igraph_bool_t), nv, perm, 
                get_igraph_vector_bool_t, set_igraph_vector_bool_t);
        SETVABV(tree, STR(a->boolvattr_names, i), &a->boolvattrs[i]);
    }

    for (i = 0; i < a->nstreattr; ++i) {
        permute(&a->streattrs[i], BUFSIZ, nv, perm, get_igraph_strvector_t,
                set_igraph_strvector_t);
        for (j = 0; j < nv; ++j) {
            if (edge[j] != -1) {
                SETEAS(tree, STR(a->streattr_names, i), edge[j], 
                       STR(a->streattrs[i], j));
            }
        }
    }

    for (i = 0; i < a->nnumeattr; ++i) {
        permute(&a->numeattrs[i], sizeof(igraph_real_t), nv, perm, 
                get_igraph_vector_t, set_igraph_vector_t);
        for (j = 0; j < nv; ++j) {
            if (edge[j] != -1) {
                SETEAN(tree, STR(a->numeattr_names, i), edge[j], 
                       VECTOR(a->numeattrs[i])[j]);
            }
        }
    }

    for (i = 0; i < a->nbooleattr; ++i) {
        permute(&a->booleattrs[i], sizeof(igraph_bool_t), nv, perm, 
                get_igraph_vector_bool_t, set_igraph_vector_bool_t);
        for (j = 0; j < nv; ++j) {
            if (edge[j] != -1) {
                SETEAN(tree, STR(a->booleattr_names, i), edge[j], 
                       VECTOR(a->booleattrs[i])[j]);
            }
        }
    }
    free(edge);
}

void _tree_attrs_destroy(tree_attrs *a)
{
    int i;

    for (i = 0; i < a->nstrvattr; ++i) {
        igraph_strvector_destroy(&a->strvattrs[i]);
    }
    for (i = 0; i < a->nnumvattr; ++i) {
        igraph_vector_destroy(&a->numvattrs[i]);
    }
    for (i = 0; i < a->nboolvattr; ++i) {
        igraph_vector_bool_destroy(&a->boolvattrs[i]);
    }

    for (i = 0; i < a->nstreattr; ++i) {
        igraph_strvector_destroy(&a->streattrs[i]);
    }
    for (i = 0; i < a->nnumeattr; ++i) {
        igraph_vector_destroy(&a->numeattrs[i]);
    }
    for (i = 0; i < a->nbooleattr; ++i) {
        igraph_vector_bool_destroy(&a->booleattrs[i]);
    }

    free(a->strvattrs);
    free(a->numvattrs);
    free(a->boolvattrs);

    free(a->streattrs);
    free(a->numeattrs);
    free(a->booleattrs);

    igraph_strvector_destroy(&a->strvattr_names);
    igraph_strvector_destroy(&a->numvattr_names);
    igraph_strvector_destroy(&a->boolvattr_names);

    igraph_strvector_destroy(&a->streattr_names);
    igraph_strvector_destroy(&a->numeattr_names);
    igraph_strvector_destroy(&a->booleattr_names);

    free(a);
}

double height(const igraph_t *tree)
{
    igraph_vector_t work;
    double ht;

    if (!igraph_ecount(tree)) {
        return 0;
    }

    igraph_vector_init(&work, 2);
    ht = _height(tree, &work, root(tree));
    igraph_vector_destroy(&work);
    return ht;
}

void write_tree_newick(const igraph_t *tree, FILE *f)
{
    // TODO: should get a more accurate estimate of the space we need
    char *s = calloc(igraph_vcount(tree) * 100, sizeof(char));
    igraph_vector_t work;

    igraph_vector_init(&work, 0);
    _write_tree_newick(tree, s, root(tree), &work);
    fprintf(f, "%s;", s);

    igraph_vector_destroy(&work);
    free(s);
}

double scale_branches(igraph_t *tree, scaling mode)
{
    int i;
    double scale;
    igraph_vector_t vec;
    if (!igraph_ecount(tree)) {
        return 1;
    }

    igraph_vector_init(&vec, igraph_vcount(tree));

    EANV(tree, "length", &vec);
    switch (mode)
    {
        case MEAN:
            scale = gsl_stats_mean(VECTOR(vec), 1, igraph_ecount(tree));
            break;
        case MEDIAN:
            igraph_vector_sort(&vec);
            scale = gsl_stats_median_from_sorted_data(VECTOR(vec), 1, igraph_ecount(tree));
            break;
        case MAX:
            scale = igraph_vector_max(&vec);
            break;
        default:
            scale = 1.;
            break;
    }

    for (i = 0; i < igraph_ecount(tree); ++i) {
        SETEAN(tree, "length", i, EAN(tree, "length", i) / scale);
    }
    igraph_vector_destroy(&vec);
    return 1.0 / scale;
}

void cut_at_time(igraph_t *tree, double t, int extant_only)
{
    int i;
    igraph_vector_t work, to_delete;
    igraph_vector_init(&work, 0);
    igraph_vector_init(&to_delete, 0);

    _cut_at_time(tree, t, root(tree), 0., extant_only, &work, &to_delete);
    igraph_vector_sort(&to_delete);
    igraph_delete_vertices(tree, igraph_vss_vector(&to_delete));
    collapse_singles(tree);

    igraph_vector_destroy(&work);
    igraph_vector_destroy(&to_delete);
}

void collapse_singles(igraph_t *tree)
{
    int i, edge;
    double bl = 0.;
    igraph_vector_t vdel, eadd, work, branch_length;

    igraph_vector_init(&vdel, 0);
    igraph_vector_init(&eadd, 0);
    igraph_vector_init(&work, 0);
    igraph_vector_init(&branch_length, 0);

    _collapse_singles(tree, root(tree), &vdel, &eadd, &branch_length, &work, &bl);
    igraph_add_edges(tree, &eadd, 0);
    for (i = 0; i < igraph_vector_size(&eadd)/2; ++i)
    {
        igraph_get_eid(tree, &edge, VECTOR(eadd)[2*i], VECTOR(eadd)[2*i+1], 1, 1);
        SETEAN(tree, "length", edge, VECTOR(branch_length)[i]);
    }
    igraph_delete_vertices(tree, igraph_vss_vector(&vdel));

    igraph_vector_destroy(&branch_length);
    igraph_vector_destroy(&vdel);
    igraph_vector_destroy(&eadd);
    igraph_vector_destroy(&work);
}

void subsample_tips(igraph_t *tree, int ntip, const gsl_rng *rng)
{
    int i, j, orig_ntip = (igraph_vcount(tree) + 1)/2;
    igraph_vector_t tips, keep_tips, keep_all, drop, degree;
    igraph_vector_ptr_t nbhd;
    igraph_vector_t *elem;

    if (ntip <= 0 || orig_ntip <= ntip) {
        return;
    }

    igraph_vector_init(&tips, 0);
    igraph_vector_init(&keep_tips, ntip);
    igraph_vector_init(&keep_all, 0);
    igraph_vector_init(&degree, 0);
    igraph_vector_init(&drop, 0);
    igraph_vector_ptr_init(&nbhd, 0);

    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        if ((int) VECTOR(degree)[i] == 0) {
            igraph_vector_push_back(&tips, i);
        }
    }
    gsl_ran_choose(rng, VECTOR(keep_tips), ntip, VECTOR(tips), orig_ntip, sizeof(igraph_real_t));

    igraph_neighborhood(tree, &nbhd, igraph_vss_vector(&keep_tips), INT_MAX, IGRAPH_IN, 0);

    for (i = 0; i < igraph_vector_ptr_size(&nbhd); ++i) {
        elem = (igraph_vector_t *) igraph_vector_ptr_e(&nbhd, i);
        for (j = 0; j < igraph_vector_size(elem); ++j) {
            igraph_vector_push_back(&keep_all, VECTOR(*elem)[j]);
        }
        igraph_vector_destroy(elem);
    }

    igraph_vector_sort(&keep_all);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        if (!igraph_vector_binsearch2(&keep_all, i))
            igraph_vector_push_back(&drop, i);
    }

    igraph_delete_vertices(tree, igraph_vss_vector(&drop));
    collapse_singles(tree);

    igraph_vector_destroy(&tips);
    igraph_vector_destroy(&keep_tips);
    igraph_vector_destroy(&keep_all);
    igraph_vector_destroy(&degree);
    igraph_vector_destroy(&drop);
    igraph_vector_ptr_destroy_all(&nbhd);
}

void subsample_tips_peerdriven(igraph_t *tree, const igraph_t *net, double p, 
        double a, int ntip, const gsl_rng *rng)
{
    int i, j, s, tip, nt = igraph_vcount(tree);
    int *tip_map = malloc(igraph_vcount(net) * sizeof(int));
    int *node_map = malloc(nt * sizeof(int));
    int *sampled = malloc(nt * sizeof(int));
    double *prob = malloc(nt * sizeof(double));
    int *idx = malloc(nt * sizeof(int));
    igraph_adjlist_t al;
    igraph_vector_t degree, keep_tips, keep_all, drop, *elem;
    igraph_vector_int_t *peers;
    igraph_vector_ptr_t nbhd;

    if (ntip <= 0 || (nt + 1) / 2 <= ntip) {
        return;
    }

    igraph_vector_init(&degree, 0);
    igraph_vector_init(&keep_tips, 0);
    igraph_vector_init(&keep_all, 0);
    igraph_vector_init(&drop, 0);
    igraph_vector_ptr_init(&nbhd, 0);
    igraph_adjlist_init(net, &al, IGRAPH_ALL);

    igraph_degree(net, &degree, igraph_vss_all(), IGRAPH_ALL, 0);
    _make_maps(tree, net, tip_map, node_map);
    for (i = 0; i < nt; ++i) {
        sampled[i] = (node_map[i] == -1 ? 1 : 0);
        idx[i] = i;
    }

    for (s = 0; s < ntip; ++s) {
        memset(prob, 0, nt * sizeof(double));
        for (i = 0; i < nt; ++i) {
            if (!sampled[i]) {
                prob[i] = p;
                peers = igraph_adjlist_get(&al, node_map[i]);
                for (j = 0; j < VECTOR(degree)[node_map[i]]; ++j) {
                    if (tip_map[VECTOR(*peers)[j]] != -1 && sampled[tip_map[VECTOR(*peers)[j]]]) {
                        prob[i] = p + a;
                        break;
                    }
                }
            }
        }
        do {
            sample_weighted(idx, &tip, 1, nt, sizeof(int), prob, 0, rng);
        } while (node_map[tip] == -1);
        sampled[tip] = 1;
        igraph_vector_push_back(&keep_tips, tip);
    }

    igraph_neighborhood(tree, &nbhd, igraph_vss_vector(&keep_tips), INT_MAX, IGRAPH_IN, 0);

    for (i = 0; i < igraph_vector_ptr_size(&nbhd); ++i) {
        elem = (igraph_vector_t *) igraph_vector_ptr_e(&nbhd, i);
        for (j = 0; j < igraph_vector_size(elem); ++j) {
            igraph_vector_push_back(&keep_all, VECTOR(*elem)[j]);
        }
        igraph_vector_destroy(elem);
    }

    igraph_vector_sort(&keep_all);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        if (!igraph_vector_binsearch2(&keep_all, i))
            igraph_vector_push_back(&drop, i);
    }

    igraph_delete_vertices(tree, igraph_vss_vector(&drop));
    collapse_singles(tree);

    igraph_vector_destroy(&drop);
    igraph_vector_destroy(&keep_tips);
    igraph_vector_destroy(&keep_all);
    igraph_vector_destroy(&degree);
    igraph_vector_ptr_destroy_all(&nbhd);
    igraph_adjlist_destroy(&al);
    free(tip_map);
    free(node_map);
    free(sampled);
    free(prob);
    free(idx);
}

void subsample(igraph_t *tree, int ntime, const double *prop, const double *t, gsl_rng *rng)
{
    int i, j, sample, v, nextant;
    int *sample_order = malloc(ntime * sizeof(int));
    double *dpth = NULL;
    igraph_vector_ptr_t nbhd;
    igraph_vector_t extant, drop, drop_all, *elem;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    igraph_vector_int_t *adj, *inc;

    igraph_vector_init(&extant, 0);
    igraph_vector_init(&drop, igraph_vcount(tree));
    igraph_vector_init(&drop_all, 0);
    igraph_vector_ptr_init(&nbhd, 0);

    // order the sampling times chronologically
    order(t, sample_order, sizeof(double), ntime, compare_doubles);

    // go through each sampling time, choosing the required number of tips
    for (i = 0; i <= ntime; ++i) {
        igraph_adjlist_init(tree, &adjlist, IGRAPH_IN);
        igraph_inclist_init(tree, &inclist, IGRAPH_IN);

        // get the depths of each node
        dpth = safe_realloc(dpth, igraph_vcount(tree) * sizeof(double));
        depths(tree, 1, dpth);
        
        sample = sample_order[i == ntime ? i - 1 : i];
        igraph_vector_clear(&extant);
        igraph_vector_clear(&drop);
        igraph_vector_clear(&drop_all);
        nextant = 0;

        // find all the nodes which are extant at the sampling time
        // by extant we mean that they are later than the sampling time, but
        // their parent is earlier
        for (v = 0; v < igraph_vcount(tree); ++v) {
            adj = igraph_adjlist_get(&adjlist, v);

            // skip the root
            if (igraph_vector_int_size(adj) == 0) {
                continue;
            }

            if (dpth[v] >= t[sample] && dpth[VECTOR(*adj)[0]] <= t[sample]) {
                igraph_vector_push_back(&extant, v);
                ++nextant;
            }
        }

        if (i < ntime) 
        {
            // choose the required number of nodes to be sampled
            igraph_vector_resize(&drop, (int) (prop[sample] * nextant));
            gsl_ran_choose(rng, VECTOR(drop), (int) (prop[sample] * nextant), 
                           VECTOR(extant), nextant, sizeof(igraph_real_t));
    
            // adjust the branch lengths of the sampled nodes
            for (j = 0; j < igraph_vector_size(&drop); ++j) {
                v = VECTOR(drop)[j];
                inc = igraph_inclist_get(&inclist, v);
                SETEAN(tree, "length", VECTOR(*inc)[0], 
                       EAN(tree, "length", VECTOR(*inc)[0]) - dpth[v] + t[sample]);
            }
        }

        // at the end, delete everybody who wasn't sampled
        else 
        {
            igraph_vector_resize(&drop, nextant);
            memcpy(VECTOR(drop), VECTOR(extant), nextant * sizeof(igraph_real_t));
        }

        // find all children of the sampled nodes and delete them
        igraph_neighborhood(tree, &nbhd, igraph_vss_vector(&drop), INT_MAX, IGRAPH_OUT, 0);

        for (j = 0; j < igraph_vector_ptr_size(&nbhd); ++j) {
            elem = (igraph_vector_t *) igraph_vector_ptr_e(&nbhd, j);
            for (v = 0; v < igraph_vector_size(elem); ++v) {
                if (i == ntime || !igraph_vector_binsearch2(&drop, VECTOR(*elem)[v])) {
                    igraph_vector_push_back(&drop_all, VECTOR(*elem)[v]);
                }
            }
            igraph_vector_destroy(elem);
        }

        igraph_delete_vertices(tree, igraph_vss_vector(&drop_all));
        igraph_adjlist_destroy(&adjlist);
        igraph_inclist_destroy(&inclist);
    }

    collapse_singles(tree);

    free(dpth);
    free(sample_order);
    igraph_vector_destroy(&extant);
    igraph_vector_destroy(&drop);
    igraph_vector_destroy(&drop_all);
    igraph_vector_ptr_destroy_all(&nbhd);
}

void depths(const igraph_t *tree, int use_branch_lengths, double *depths)
{
    igraph_vector_t work;
    igraph_vector_init(&work, 1);
    _depths(tree, depths, &work, root(tree), 0.0, use_branch_lengths);
    igraph_vector_destroy(&work);
}

/* Private */

void _depths(const igraph_t *tree, double *depths, igraph_vector_t *work, 
        int root, double parent_depth, int use_branch_lengths)

{
    int lc, rc;

    igraph_incident(tree, work, root, IGRAPH_IN);
    if (igraph_vector_size(work) > 0) {
        if (use_branch_lengths) {
            depths[root] = parent_depth + EAN(tree, "length", (int) VECTOR(*work)[0]);
        }
        else {
            depths[root] = parent_depth + 1;
        }
    }
    else {
        depths[root] = parent_depth;
    }

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) > 0)
    {
        lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
        _depths(tree, depths, work, lc, depths[root], use_branch_lengths);
        _depths(tree, depths, work, rc, depths[root], use_branch_lengths);
    }
}

int _ladderize(igraph_t *tree, igraph_vector_t *work, int root, int *perm)
{
    igraph_es_t es;
    int do_swap, lc, rc, lsize, rsize;

    // if a leaf, just add the next number
    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        perm[0] = root;
        return 1;
    }

    // get the children
    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    lc = (int) VECTOR(*work)[0];
    rc = (int) VECTOR(*work)[1];

    // order first by size
    lsize = _ladderize(tree, work, lc, perm);
    rsize = _ladderize(tree, work, rc, &perm[lsize]);
    do_swap = (lsize > rsize);

    // order next by branch length
    if (lsize == rsize)
    {
        igraph_es_incident(&es, root, IGRAPH_OUT);
        igraph_cattribute_EANV(tree, "length", es, work);
        do_swap = (VECTOR(*work)[0] > VECTOR(*work)[1]);
    }

    // swap order of children if necessary
    if (do_swap) {
        rotl(perm, (rsize + lsize) * sizeof(int), lsize * sizeof(int));
    }

    // now add the root
    perm[lsize + rsize] = root;
    return lsize + rsize + 1;
}

double _height(const igraph_t *tree, igraph_vector_t *work, int root)
{
    int lc, rc;
    double height;
    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        igraph_incident(tree, work, root, IGRAPH_IN);
        return EAN(tree, "length", (int) VECTOR(*work)[0]);
    }

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
    height = fmax(_height(tree, work, lc), _height(tree, work, rc));

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_IN, 0);
    if ((int) VECTOR(*work)[0] == 0)
        return height;
    igraph_incident(tree, work, root, IGRAPH_IN);
    return EAN(tree, "length", (int) VECTOR(*work)[0]) + height;
}

int _write_tree_newick(const igraph_t *tree, char *out, int root,
        igraph_vector_t *work)
{
    double length;
    int nchar, left_child, right_child, is_root = 0;

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_IN, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        length = 0.;
        is_root = 1;
    }
    else
    {
        igraph_incident(tree, work, root, IGRAPH_IN);
        length = EAN(tree, "length", VECTOR(*work)[0]);
    }

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        if (igraph_cattribute_has_attr(tree, IGRAPH_ATTRIBUTE_VERTEX, "id"))
            return sprintf(out, "%s:%f", VAS(tree, "id", root), length);
        else
            return sprintf(out, "%d:%f", root, length);
    }

    igraph_neighbors(tree, work, root, IGRAPH_OUT);

    left_child = VECTOR(*work)[0];
    right_child = VECTOR(*work)[1];

    nchar = sprintf(out, "(");
    nchar += _write_tree_newick(tree, &out[nchar], left_child, work);
    nchar += sprintf(&out[nchar], ",");
    nchar += _write_tree_newick(tree, &out[nchar], right_child, work);
    nchar += sprintf(&out[nchar], ")");

    if (igraph_cattribute_has_attr(tree, IGRAPH_ATTRIBUTE_VERTEX, "id"))
        nchar += sprintf(&out[nchar], "%s:%f", VAS(tree, "id", root), length);
    else
        nchar += sprintf(&out[nchar], "%d:%f", root, length);
    return nchar;
}

void _cut_at_time(igraph_t *tree, double t, int root, double troot, 
        int extant_only, igraph_vector_t *work, igraph_vector_t *to_delete)
{
    int i, lc, rc;
    double tnode;

    igraph_incident(tree, work, root, IGRAPH_IN);
    if (igraph_vector_size(work) > 0) {
        tnode = EAN(tree, "length", VECTOR(*work)[0]);
    }
    else {
        tnode = 0.;
    }

    if (tnode + troot > t)
    {
        // adjust my branch length
        SETEAN(tree, "length", VECTOR(*work)[0], t - troot);

        // delete all descendents
        igraph_subcomponent(tree, work, root, IGRAPH_OUT);
        for (i = 0; i < igraph_vector_size(work); ++i) 
        {
            if ((int) VECTOR(*work)[i] != root) {
                igraph_vector_push_back(to_delete, VECTOR(*work)[i]);
            }
        }
    }
    else
    {
        igraph_neighbors(tree, work, root, IGRAPH_OUT);
        if (igraph_vector_size(work) > 0)
        {
            lc = (int) VECTOR(*work)[0];
            rc = (int) VECTOR(*work)[1];

            _cut_at_time(tree, t, lc, troot + tnode, extant_only, work, to_delete);
            _cut_at_time(tree, t, rc, troot + tnode, extant_only, work, to_delete);
        }

        // if we're sampling extant nodes at time t only, and this node isn't
        // extant at time t, then we delete it
        else if (extant_only && troot + tnode < t)
        {
            igraph_vector_push_back(to_delete, (igraph_real_t) root);
        }
    }
}

int _collapse_singles(igraph_t *tree, int root, igraph_vector_t *vdel,
        igraph_vector_t *eadd, igraph_vector_t *branch_length, 
        igraph_vector_t *work, double *bl)
{
    int lc, rc, new_lc, new_rc;
    double lbl, rbl;

    igraph_incident(tree, work, root, IGRAPH_IN);
    if (igraph_vector_size(work) > 0)
        bl[0] = EAN(tree, "length", (int) VECTOR(*work)[0]);
    else
        bl[0] = 0;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) == 0)
    {
        return root;
    }
    else if (igraph_vector_size(work) == 1)
    {
        igraph_vector_push_back(vdel, root);
        lc = (int) VECTOR(*work)[0];
        new_lc = _collapse_singles(tree, lc, vdel, eadd, branch_length, work, &lbl);
        bl[0] += lbl;
        return new_lc;
    }
    else
    {
        lc = (int) VECTOR(*work)[0];
        rc = (int) VECTOR(*work)[1];
        new_lc = _collapse_singles(tree, lc, vdel, eadd, branch_length, work, &lbl);
        new_rc = _collapse_singles(tree, rc, vdel, eadd, branch_length, work, &rbl);

        if (new_lc != lc)
        {
            igraph_vector_push_back(eadd, root);
            igraph_vector_push_back(eadd, new_lc);
            igraph_vector_push_back(branch_length, lbl);
        }
        if (new_rc != rc)
        {
            igraph_vector_push_back(eadd, root);
            igraph_vector_push_back(eadd, new_rc);
            igraph_vector_push_back(branch_length, rbl);
        }

        return root;
    }
}

void _get_node_ids(const igraph_t *g, igraph_strvector_t *ids)
{
    igraph_attribute_type_t id_type = get_igraph_id_type(g);
    char buf[BUFSIZ];
    int i;

    for (i = 0; i < igraph_vcount(g); ++i) {
        switch (id_type) {
            case IGRAPH_ATTRIBUTE_NUMERIC:
                sprintf(buf, "%d", (int) VAN(g, "id", i));
                break;
            case IGRAPH_ATTRIBUTE_STRING:
                sprintf(buf, "%s", VAS(g, "id", i));
                break;
            default:
                sprintf(buf, "%d", i);
                break;
        }
        igraph_strvector_add(ids, buf);
    }
}

void _make_maps(const igraph_t *tree, const igraph_t *net, int *tip_map, int *node_map)
{
    int i;
    igraph_vector_t degree;
    igraph_strvector_t tree_ids, net_ids;

    igraph_vector_init(&degree, 0);
    igraph_strvector_init(&tree_ids, 0);
    igraph_strvector_init(&net_ids, 0);

    _get_node_ids(tree, &tree_ids);
    _get_node_ids(net, &net_ids);

    // set all the non-tip nodes' ids to null
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        if ((int) VECTOR(degree)[i] > 0) {
            igraph_strvector_set(&tree_ids, i, "\0");
        }
    }
    match(&tree_ids, &net_ids, node_map, BUFSIZ, igraph_vcount(tree),
          igraph_vcount(net), get_igraph_strvector_t, compare_strings);
    match(&net_ids, &tree_ids, tip_map, BUFSIZ, igraph_vcount(net),
          igraph_vcount(tree), get_igraph_strvector_t, compare_strings);

    igraph_vector_destroy(&degree);
    igraph_strvector_destroy(&tree_ids);
    igraph_strvector_destroy(&net_ids);
}
