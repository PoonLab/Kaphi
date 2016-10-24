/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>  // provides NEW_INTEGER and NEW_NUMERIC

#include "rinterface.h"
#include "treekernel/tree.h"

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
    SEXP myint;
    int * p_myint;
    int len = 1;

    R_SEXP_to_igraph(graph, &g);
    PROTECT(myint = NEW_INTEGER(len));
    p_myint = INTEGER_POINTER(myint);
    p_myint[0] = (int)igraph_vcount(&g);
    igraph_destroy(&g);  // release allocated memory

    UNPROTECT(1);
    return myint;
}

SEXP R_Kaphi_get_edge_lengths(SEXP graph) {
    /* check that graph has edge length attributes */
    igraph_t g;
    igraph_vector_t edge_lengths;
    SEXP res;
    long int len = (int) igraph_ecount(&g);

    R_SEXP_to_igraph(graph, &g);
    igraph_vector_init(&edge_lengths, len);

    PROTECT(res=NEW_NUMERIC(igraph_vector_size(&edge_lengths)));

    igraph_bool_t has_edge_lengths = igraph_cattribute_has_attr(&g, IGRAPH_ATTRIBUTE_EDGE, "length");
    if (has_edge_lengths) {
        EANV(&g, "length", &edge_lengths);
        for (int i = 0; i < len; i++) {
            REAL(res)[i] = igraph_vector_e(&edge_lengths, i);
        }
    }
    igraph_vector_destroy(&edge_lengths);
    igraph_destroy(&g);
    UNPROTECT(1);
    return res;
}
