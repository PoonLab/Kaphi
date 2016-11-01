/* modeled on rigraph/tools/stimulus/rinterface_extra.c */

#include <igraph.h>
#include <igraph_error.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>  // provides NEW_INTEGER and NEW_NUMERIC

#include "rinterface.h"
#include "treekernel/tree.h"
#include <assert.h>
#include <Judy.h>

// function prototypes
int *production(const igraph_t *tree);
int *children(const igraph_t *tree);
double *branch_lengths(const igraph_t *tree);





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
    /* check that graph from R has edge length attributes */
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


void import_R_igraph(SEXP graph, igraph_t *g, igraph_vector_t * edge_lengths) {
    igraph_es_t edge_selector;
    igraph_eit_t edge_iterator;
    int len;
    double * temp;

    R_SEXP_to_igraph(graph, g);
    len = igraph_ecount(g);
    igraph_vector_init(edge_lengths, len);

    if (R_igraph_attribute_has_attr(g, IGRAPH_ATTRIBUTE_EDGE, "length")) {
        // retrieve attribute for all edges in order of edge ID's
        igraph_es_all(&edge_selector, IGRAPH_EDGEORDER_ID);
        R_igraph_attribute_get_numeric_edge_attr(g, "length",
            edge_selector, edge_lengths
        );
        // FIXME: can't get SETEAN to work here
    }
}


SEXP R_Kaphi_kernel(SEXP graph1, SEXP graph2, SEXP arg_lambda, SEXP arg_sigma, SEXP arg_rho) {
    // unpack real-valued arguments
    double lam = REAL(arg_lambda)[0];
    double sig = REAL(arg_sigma)[0];
    double rho = REAL(arg_rho)[0];
    int i, c1, c2, n1, n2, coord, npairs;
    int * production1, * production2, * children1, * children2;
    double val, tmp, K = 0;
    SEXP result;
    igraph_es_t edge_selector;
    igraph_vector_t e1, e2;

    Pvoid_t delta = (Pvoid_t) NULL;  // Judy1 array
    Word_t bytes = 0;

    igraph_t t1, t2;

    assert(decay_factor > 0.0 && decay_factor <= 1.0);
    assert(rbf_variance > 0.0);
    assert(igraph_vcount(t1) < 65535);

    import_R_igraph(graph1, &t1, &e1);
    import_R_igraph(graph2, &t1, &e2);

    // extract production states of internal nodes in the tree
    production1 = production(&t1);
    production2 = production(&t2);

    children1 = children(&t1);
    children2 = children(&t2);


    PROTECT(result = NEW_NUMERIC(1));
    REAL(result)[0] = 9.0;
    UNPROTECT(1);

    return result;
}


SEXP R_Kaphi_get_productions(SEXP graph) {
    // for unit test of production()
    igraph_t g;
    int * p;
    SEXP result;
    int nnode;
    int * p_result;

    R_SEXP_to_igraph(graph, &g);
    nnode = igraph_vcount(&g);

    PROTECT(result = NEW_INTEGER(nnode));
    p_result = INTEGER_POINTER(result);

    p = production(&g);
    for (int i = 0; i < nnode; i++) {
        p_result[i] = p[i];
    }
    UNPROTECT(1);
    return result;
}


int * production(const igraph_t * tree)
{
    /*
        Get productions for each node.  Productions take the following values:
        0 = node is terminal
        1 = node emits two internal nodes
        2 = node emits one terminal node
        3 = node emits two terminal nodes
    */
    int i, nnode = igraph_vcount(tree);
    int *p = malloc(nnode * sizeof(int));
    igraph_vector_t vec, nbr;
    igraph_vector_init(&nbr, 2);
    igraph_vector_init(&vec, 2);

    // retrieves out-degree size of each vertex
    // igraph_vss_all is a selector for all vertices in increasing order of vertex ID
    // vertex IDs are assigned when phylo object is converted to igraph in R
    igraph_degree(tree, &vec, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < nnode; ++i) {
        if ((int) VECTOR(vec)[i] == 0) {
            // terminal node has no productions
            // TODO: this is where we have to deal with labels
            p[i] = 0;
        } else {
            // retrieves IDs of all out-adjacent nodes of i-th node
            igraph_neighbors(tree, &nbr, i, IGRAPH_OUT);

            // use neighbor IDs to index into out-degree vector
            // if result is 0, then neighbor is terminal
            p[i] = ((int) VECTOR(vec)[(int) VECTOR(nbr)[0]] == 0) +
                   ((int) VECTOR(vec)[(int) VECTOR(nbr)[1]] == 0) + 1;
        }
    }
    igraph_vector_destroy(&nbr);
    igraph_vector_destroy(&vec);
    return p;
}


SEXP R_Kaphi_get_children(SEXP graph) {
    // for unit test of production()
    igraph_t g;
    int * c;
    SEXP result;
    int nnode;
    int * p_result;

    R_SEXP_to_igraph(graph, &g);

    // allocate double length to store two child indices per node
    nnode = igraph_vcount(&g) * 2;

    PROTECT(result = NEW_INTEGER(nnode));
    p_result = INTEGER_POINTER(result);

    c = children(&g);
    for (int i = 0; i < nnode; i++) {
        p_result[i] = c[i];
    }
    UNPROTECT(1);
    return result;
}

int *children(const igraph_t *tree)
{
    /* get indices of children from each node */
    igraph_adjlist_t al;
    igraph_vector_int_t *nbr;
    int i;
    int *children = malloc(igraph_vcount(tree) * 2 * sizeof(int));

    // constructs out-edge adjacency list from graph
    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        // retrieve out-edge adjacency (child node) list for i-th node
        nbr = igraph_adjlist_get(&al, i);
        if (igraph_vector_int_size(nbr) > 0) {
            // we assume the tree is binary, only ever two children
            children[2*i] = VECTOR(*nbr)[0];
            children[2*i+1] = VECTOR(*nbr)[1];
        } else {
            // store negative index to indicate no children
            children[2*i] = children[2*i+1] = -1;
        }
    }

    igraph_adjlist_destroy(&al);
    return children;
}


SEXP R_Kaphi_get_branch_lengths(SEXP graph) {
    // for unit test of production()
    igraph_t g;
    double * bl;
    SEXP result;
    int nnode;
    igraph_vector_t edge_lengths;

    // turn on attribute handling for C igraphs
    // igraph_i_set_attribute_table(&igraph_cattribute_table);

    import_R_igraph(graph, &g, &edge_lengths);
    for (int i=0; i<igraph_ecount(&g); i++) {
        fprintf(stdout, "%f\n", igraph_vector_e(&edge_lengths, i));
    }

    // allocate double length to store two child indices per node
    nnode = igraph_vcount(&g) * 2;

    PROTECT(result = NEW_NUMERIC(nnode));

    //bl = branch_lengths(&g);
    for (int i = 0; i < nnode; i++) {
        //REAL(result)[i] = bl[i];
    }
    UNPROTECT(1);
    return result;
}

