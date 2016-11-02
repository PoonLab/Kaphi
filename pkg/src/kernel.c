/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>
#include <igraph_error.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>  // provides NEW_INTEGER and NEW_NUMERIC

//#include "rinterface.h"
#include "newick_parser.h"
#include "tree.h"


igraph_t * R_Kaphi_parse_newick(SEXP newick) {
    // derived from Rosemary's parse_newick() function
    igraph_t *tree;
    const char * newick_str = CHAR(STRING_ELT(newick, 0));

    extern int yynode;  // in newick_parser.y
    igraph_vector_t edge, branch_length, size;
    igraph_strvector_t label;

    // initialize vector containers
    igraph_vector_init(&edge, 0);
    igraph_vector_init(&size, 0);
    igraph_vector_init(&branch_length, 0);
    igraph_strvector_init(&label, 0);

    yynode = 0;
    yy_scan_string(newick_str);
    yyparse(&edge, &size, &branch_length, &label);  // attempt to parse input stream

    tree = malloc(sizeof(igraph_t));
    igraph_empty(tree, igraph_vector_size(&size), 1);
    igraph_add_edges(tree, &edge, 0);

    for (int i = 0; i < igraph_vector_size(&size); ++i)
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


