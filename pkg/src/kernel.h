typedef struct {
    igraph_strvector_t *strvattrs;
    igraph_vector_t *numvattrs;
    igraph_vector_bool_t *boolvattrs;

    igraph_strvector_t *streattrs;
    igraph_vector_t *numeattrs;
    igraph_vector_bool_t *booleattrs;

    int nstrvattr;
    int nnumvattr;
    int nboolvattr;

    int nstreattr;
    int nnumeattr;
    int nbooleattr;

    igraph_strvector_t strvattr_names;
    igraph_strvector_t numvattr_names;
    igraph_strvector_t boolvattr_names;

    igraph_strvector_t streattr_names;
    igraph_strvector_t numeattr_names;
    igraph_strvector_t booleattr_names;
} tree_attrs;

tree_attrs *_get_tree_attrs(const igraph_t *tree);

typedef enum {
    MEAN,
    MEDIAN,
    MAX,
    NONE
} scaling;
