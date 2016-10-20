/** \file util.h
 * \brief Utility functions. 
 */

#ifndef UTIL_H
#define UTIL_H

#include <gsl/gsl_rng.h>

#include "../igraph/include/igraph.h"

/** Set seeds for all random number generators.
 *
 * This sets the seeds of the random number generators used by libc, GSL, and
 * igraph all to the same value. Since GSL's random generation requires that
 * you pass the generator as a parameter, the generator is returned.
 *
 * \param[in] seed the new seed
 * \return the GSL random generator
 */
gsl_rng *set_seed(int seed);

/** Compare two doubles.
 *
 * This is intended for use as a comparator in qsort and order.
 * 
 * \param[in] a,b doubles to compare
 * \return an integer with the same sign as a - b, or zero if a == b
 * \sa order()
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
int compare_doubles (const void * a, const void * b);

/** Compare two integers.
 *
 * This is intended for use as a comparator in qsort and order.
 * 
 * \param[in] a,b integers to compare
 * \return an integer with the same sign as a - b, or zero if a == b
 * \sa order()
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
int compare_ints (const void * a, const void * b);

/** Compare two strings.
 *
 * This is intended for use as a comparator in qsort and order. It's a thin
 * wrapper around strcmp.
 * 
 * \param[in] a,b strings to compare
 * \return 0 if a == b, non-zero otherwise
 * \sa strcmp()
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
int compare_strings(const void * a, const void *b);

/** Get the sorting order of an array.
 *
 * This function returns a permutation of x which would put x in sorted order.
 * It's like the order function in R, only zero-indexed. The compar function is
 * like one which is used for qsort.
 *
 * \param[in] base an array of items to order
 * \param[out] order the result will be stored here
 * \param[in] nitems number of items to order
 * \param[in] compar function which compares two elements of base
 *
 * \sa https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
void order(const void *base, int *order, size_t size, int nitems,
        int (*compar) (const void *, const void *));

/** Circular left shift a block of memory.
 *
 * \param[in] x the block of memory to shift
 * \param[in] nx the size of x, in bytes
 * \param[in] n the number of bytes to shift left
 */
void rotl(void *x, size_t nx, size_t n);

/** Get a log10 scale factor for some numbers.
 *
 * The numbers can be divided by 10 to the power of the returned value, to keep
 * them in a reasonable range.
 *
 * \param[in] x the numbers to scale
 * \param[in] n the length of x
 * \return the base 10 logarithm of a scaling factor to divide both numbers by
 */
int get_scale(double *x, int n);

/** Find the index of the maximum element in an array of doubles.
 *
 * \param[in] x an array of doubles
 * \param[in] n the length of x
 * \return the index of the largest element in x
 */
int which_max(double *x, int n);

/** Sum the elements of a double array.
 *
 * \param[in] x an array of doubles
 * \param[in] n the number of elements to sum, up to the length of x
 * \return the sum of the first n elements of x
 */
double sum_doubles(const double *x, int n);

/** Find the maximum element in an array of doubles.
 *
 * \param[in] x an array of doubles
 * \param[in] n the length of x
 * \return the largest element in x
 */
double max_doubles(const double *x, int n);

/** Allocate a block of memory, or abort if out of memory.
 *
 * \param[in,out] ptr place to put new memory block
 * \param[in] size amount of memory to allocate
 */
void *safe_realloc(void *ptr, size_t size);

/** Permute an array in place
 *
 * You have to supply a getter and setter for the array. This seems a bit
 * unnecessary when using regular arrays, but it's so that the function can
 * also be applied to eg. igraph_strvector_t which can't be indexed with []
 * syntax.
 *
 * \param[in,out] v array to permute
 * \param[in] size size of array elements
 * \param[in] nitems number of elements in array
 * \param[in] perm permutation to apply
 * \param[in] get returns a pointer to the nth element of v
 * \param[in] set sets the nth element of v
 */
void permute(void *v, size_t size, int nitems, const int *perm, 
             void (*get) (const void *, int, void *),
             void (*set) (void *, int, const void *));

/** Match the elements of two arrays.
 *
 * This is like the match function in R. Any values in x which don't have a
 * match in table will have their corresponding entry in pos set to -1. Like in
 * permute, you have to supply a getter for the array.
 *
 * \param[in] x the values to be matched
 * \param[in] table the values to be matched against
 * \param[out] pos the first positions of x in table will be placed here
 * \param[in] size size of elements in x
 * \param[in] nx number of elements in x
 * \param[in] ntable number of elements in table
 * \param[in] perm permutation to apply
 * \param[in] get returns a pointer to the nth element of v
 * \param[in] compar returns zero if and only if its arguments are equal
 * \sa https://stat.ethz.ch/R-manual/R-devel/library/base/html/match.html
 */
void match(const void *x, const void *table, int *pos, size_t size, int nx, int ntable,
           void (*get) (const void *, int, void *),
           int (*compar) (const void *, const void *));

/** Take a weighted sample of array elements.
 *
 * \param[in] x elements to sample
 * \param[out] dest sampled elements will be placed here
 * \param[in] k number of elements to sample
 * \param[in] n number of elements in x
 * \param[in] size size of each element in x
 * \param[in] prob sampling weights of each element in x
 * \param[in] replace whether to sample with replacement
 * \param[in] rng GSL random number generator
 */
void sample_weighted(const void *x, void *dest, size_t k, size_t n, size_t size, 
                     const double *prob, int replace, const gsl_rng *rng);

/** Set the nth element of an igraph_vector_t
 *
 * This provides an alternative inteface to igraph_vector_set to be useable in
 * permute().
 *
 * \param[in,out] v vector to be modified
 * \param[in] n index of v to set
 * \param[in] value pointer to new value of v[n]
 */
void set_igraph_vector_t(void *v, int n, const void *value);

/** Get the nth element of an igraph_vector_t
 *
 * This provides an alternative inteface to igraph_strvector_get to be useable
 * in permute().
 *
 * \param[in,out] v vector to be queried
 * \param[in] n index of v to get
 * \param[out] value v[n] will be stored here
 */
void get_igraph_vector_t(const void *v, int n, void *value);

/** Set the nth element of an igraph_vector_bool_t
 *
 * This provides an alternative inteface to igraph_vector_bool_set to be
 * useable in permute().
 *
 * \param[in,out] v vector to be modified
 * \param[in] n index of v to set
 * \param[in] value pointer to new value of v[n]
 */
void set_igraph_vector_bool_t(void *v, int n, const void *value);

/** Get the nth element of an igraph_vector_bool_t
 *
 * This provides an alternative inteface to igraph_strvector_bool_get to be
 * useable in permute().
 *
 * \param[in,out] v vector to be queried
 * \param[in] n index of v to get
 * \param[out] value v[n] will be stored here
 */
void get_igraph_vector_bool_t(const void *v, int n, void *value);

/** Set the nth element of an igraph_strvector_t
 *
 * This provides an alternative inteface to igraph_strvector_set to be useable
 * in permute().
 *
 * \param[in,out] v vector to be modified
 * \param[in] n index of v to set
 * \param[in] value pointer to new value of v[n]
 */
void set_igraph_strvector_t(void *v, int n, const void *value);

/** Get the nth element of an igraph_strvector_t
 *
 * This provides an alternative inteface to igraph_strvector_set to be useable
 * in permute().
 *
 * \param[in,out] v vector to be queried
 * \param[in] n index of v to get
 * \param[out] value v[n] will be stored here
 */
void get_igraph_strvector_t(const void *v, int n, void *value);

/** Find the type of ID a graph has.
 * 
 * If IGRAPH_ATTRIBUTE_NUMERIC or IGRAPH_ATTRIBUTE_STRING, it means there is a
 * numeric or string attribute called "id". If it is IGRAPH_ATTRIBUTE_DEFAULT,
 * it means that there are no "id" attributes, and the vertex indices should be
 * used as the ids.
 *
 * \param[in] g graph to query
 * \return an igraph_attribute_type_t corresponding to the ID type
 */
igraph_attribute_type_t get_igraph_id_type(const igraph_t *g);
#endif
