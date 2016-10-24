/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

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
    UNPROTECT(1);
    return myint;
}

SEXP R_Kaphi_get_branch_lengths(SEXP graph) {
    igraph_t g;
    SEXP my;
    int * p_myint;
    int len = 1;

    R_SEXP_to_igraph(graph, &g);
    PROTECT(myint = NEW_INTEGER(len));
    p_myint = INTEGER_POINTER(myint);
    p_myint[0] = (int)igraph_vcount(&g);
    UNPROTECT(1);
    return myint;
}