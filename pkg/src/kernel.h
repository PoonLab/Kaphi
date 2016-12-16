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
