/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>
#include <igraph_error.h>

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
