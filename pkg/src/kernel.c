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
#include "treestats.h"

void yy_scan_string(const char *);




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

    yynode = 0;  // set iterator for newick_parser
    yy_scan_string(newick_str);  // prepare to take scanner's input from string
    yyparse(&edge, &size, &branch_length, &label);  // parse input stream

    // allocate memory for igraph object
    tree = malloc(sizeof(igraph_t));
    igraph_empty(tree, igraph_vector_size(&size), 1);  // initialize with empty graph
    igraph_add_edges(tree, &edge, 0);

    for (int i = 0; i < igraph_vector_size(&size); ++i)
    {
        igraph_incident(tree, &edge, i, IGRAPH_IN);
        if (igraph_vector_size(&edge) > 0) {
            SETEAN(tree, "length", (int) VECTOR(edge)[0], VECTOR(branch_length)[i]);
        }
        SETVAS(tree, "id", i, STR(label,i));  // assign node label to vertex
    }

    // free up memory allocated to vectors
    igraph_vector_destroy(&edge);
    igraph_vector_destroy(&size);
    igraph_vector_destroy(&branch_length);
    igraph_strvector_destroy(&label);
    return tree;
}

SEXP R_Kaphi_kernel(SEXP nwk1, SEXP nwk2, SEXP lambda, SEXP sigma, SEXP rho, SEXP use_label, SEXP gamma, SEXP normalize) {
    SEXP result;

    // unpack real-valued arguments
    double decay_factor = REAL(lambda)[0];
    double gauss_factor = REAL(sigma)[0];
    double sst_control = REAL(rho)[0];
    double label_factor = REAL(gamma)[0];
    double do_normalize = REAL(normalize)[0];
    long int *new_label1 = 0;
    long int *new_label2 = 0;

    double knum, kdenom = 1.;  // numerator and denominator

    // enable C igraph attribute handler
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // parse SEXP arguments passed from R
    igraph_t * t1 = R_Kaphi_parse_newick(nwk1);
    igraph_t * t2 = R_Kaphi_parse_newick(nwk2);

    if (INTEGER(use_label)[0]) {
        new_label1 = (long int*)malloc(sizeof(long int) * (igraph_vcount(t1)));
        new_label2 = (long int*)malloc(sizeof(long int) * (igraph_vcount(t2)));
        get_labels(t1, t2, new_label1, new_label2);
    }

    // ladderize and branch scaling is handled on R side

    if (do_normalize) {
        // see Collins and Duffey, NIPS 2001
        kdenom = sqrt(kernel(t1, t1, decay_factor, gauss_factor, sst_control, new_label1, new_label1, label_factor)) *
                 sqrt(kernel(t2, t2, decay_factor, gauss_factor, sst_control, new_label2, new_label2, label_factor));
    }
    knum = kernel(t1, t2, decay_factor, gauss_factor, sst_control, new_label1, new_label2, label_factor);

    // transfer the result to a container to pass back to R
    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = knum / kdenom;
    UNPROTECT(1);

    // free memory allocated for trees
    igraph_destroy(t1);
    igraph_destroy(t2);
    if (INTEGER(use_label)[0]) {
        free(new_label1);
        free(new_label2);
    }

    return (result);
}
