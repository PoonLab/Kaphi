/*
 *  This file is part of Kaphi.
 *
 *  Kaphi is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Kaphi is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.
 */

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

SEXP R_Kaphi_nLTT(SEXP nwk1, SEXP nwk2) {
    SEXP result;

    // enable C igraph attribute handler
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // parse SEXP arguments passed from R
    igraph_t * t1 = R_Kaphi_parse_newick(nwk1);
    igraph_t * t2 = R_Kaphi_parse_newick(nwk2);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = nLTT(t1, t2);
    UNPROTECT(1);

    igraph_destroy(t1);
    igraph_destroy(t2);
    return(result);
}

SEXP R_Kaphi_sackin(SEXP nwk, SEXP arg_useBL) {
    SEXP result;
    igraph_t * t1;
    int use_branch_lengths = INTEGER(arg_useBL)[0];

    if (use_branch_lengths) {
        igraph_i_set_attribute_table(&igraph_cattribute_table);
    }
    t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = sackin(t1, use_branch_lengths);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_colless(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = colless(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_cophenetic(SEXP nwk, SEXP arg_useBL) {
    SEXP result;
    igraph_t * t1;
    int use_branch_lengths = INTEGER(arg_useBL)[0];

    if (use_branch_lengths) {
        igraph_i_set_attribute_table(&igraph_cattribute_table);
    }
    t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = cophenetic(t1, use_branch_lengths);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_ladder_length(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = ladder_length(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_il_nodes(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = il_nodes(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_width(SEXP nwk) {
    SEXP result;

    igraph_i_set_attribute_table(&igraph_cattribute_table); // <-- inserted
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = width(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_max_delta_width(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = max_delta_width(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_cherries(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = cherries(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_prop_unbalanced(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = prop_unbalanced(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_avg_unbalance(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = avg_unbalance(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_pybus_gamma(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = pybus_gamma(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}

SEXP R_Kaphi_internal_terminal_ratio(SEXP nwk) {
    SEXP result;
    igraph_t * t1 = R_Kaphi_parse_newick(nwk);

    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = internal_terminal_ratio(t1);
    UNPROTECT(1);

    igraph_destroy(t1);
    return(result);
}
