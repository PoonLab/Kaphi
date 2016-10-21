/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "rinternals.h"
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

SEXP R_Kaphi_treekernel(SEXP graph) {
    igraph_t g;
    SEXP result, dim;
    R_SEXP_to_igraph(graph, &g);

    tree_attrs * a = _get_tree_attrs()

    PROTECT(dim=NEW_INTEGER(1));
    INTEGER(dim)[0] = 42;
    SET_DIM(result, dim);
    UNPROTECT(1);
    return result;
}
