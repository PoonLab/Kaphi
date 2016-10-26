/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>
#include <igraph_error.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>  // provides NEW_INTEGER and NEW_NUMERIC

#include "rinterface.h"
#include "kernel.h"

//#include "treekernel/tree.h"
#include <assert.h>
#include <Judy.h>


SEXP R_Kaphi_test(void) {
   SEXP myint;
   int *p_myint;
   int len = 5;
   PROTECT(myint = NEW_INTEGER(len));  // Allocating storage space
   p_myint = INTEGER_POINTER(myint);
   for (int i = 0; i < len; i++) {
        p_myint[i] = 7-i;
   }
   UNPROTECT(1);
   return myint;
}


SEXP R_Kaphi_nodecount(SEXP graph) {
    igraph_t g;
    SEXP result;

    R_SEXP_to_igraph(graph, &g);
    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = igraph_vcount(&g);

    UNPROTECT(1);
    return result;
}

SEXP R_Kaphi_get_edge_lengths(SEXP graph) {
    /* check that graph has edge length attributes */
    igraph_t g;
    igraph_vector_t edge_lengths;
    long int len;
    SEXP result;
    igraph_es_t edge_selector;

    // convert R igraph to C igraph
    R_SEXP_to_igraph(graph, &g);

    len = (int) igraph_ecount(&g);
    igraph_vector_init(&edge_lengths, len);

    // maybe we need to copy attributes over?
    if(R_igraph_attribute_has_attr(&g,
        IGRAPH_ATTRIBUTE_EDGE, "length")) {

        // retrieve attribute for all edges in order of edge ID's
        igraph_es_all(&edge_selector, IGRAPH_EDGEORDER_ID);
        //EANV(&g, "length", &edge_lengths);
        R_igraph_attribute_get_numeric_edge_attr(&g, "length", edge_selector, &edge_lengths);
    } // else function will return a vector of zeroes

    PROTECT(result = NEW_NUMERIC(len));
    for (int i = 0; i < len; i++) {
        REAL(result)[i] = igraph_vector_e(&edge_lengths, i);
    }
    //igraph_vector_destroy(&edge_lengths);

    UNPROTECT(1);
    return result;
}

SEXP R_Kaphi_test2(SEXP x) {
    // can we just pass a double?
    double xx = REAL(x)[0];
    SEXP result;
    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = xx;
    UNPROTECT(1);
    return result;
}

SEXP R_Kaphi_kernel(SEXP g1, SEXP g2, SEXP arg_lambda, SEXP arg_sigma, SEXP arg_rho) {
    // unpack real-valued arguments
    double lam = REAL(arg_lambda)[0];
    double sig = REAL(arg_sigma)[0];
    double rho = REAL(arg_rho)[0];
    int i, c1, c2, n1, n2, coord, npairs;
    int *production1, *production2, *children1, *children2;
    double val, tmp, K = 0;
    SEXP result;

    Pvoid_t delta = (Pvoid_t) NULL;  // Judy1 array
    Word_t bytes = 0;

    assert(decay_factor > 0.0 && decay_factor <= 1.0);
    assert(rbf_variance > 0.0);
    assert(igraph_vcount(t1) < 65535);

    PROTECT(result = NEW_NUMERIC(1));
    UNPROTECT(1);

    return result;
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
