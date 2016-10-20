#include "igraph.h"

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include "rinterface.h"

#include <stdio.h>

SEXP R_igraph_test(SEXP graph) {
    igraph_t g;
    SEXP result;

    R_SEXP_to_igraph(graph, &g);
    PROTECT(result=NEW_INTEGER(0));
    UNPROTECT(1);
    return result;
}
