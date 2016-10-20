/** \file smc.h
 *
 * Generic approximate approximate Bayesian computation with sequential Monte
 * Carlo.
 */


#ifndef SMC_H
#define SMC_H

#include <gsl/gsl_rng.h>
#include "stats.h"

#define MAX_DIST_PARAMS 2


/** \struct smc_config
 *  \brief Configuration parameters for ABC-SMC
 *
 *  This struct contains all of the parameters needed to run ABC-SMC. To make
 *  the implementation as generic as possible, the domain-specific parts of the
 *  algorithm are passed in as functions. The first parameter to each function
 *  is a GSL random number generator, which *must* be used instead of rand() to
 *  generate random numbers. This is to ensure thread-safety, since each thread
 *  will be passed its own RNG. The GSL library has a variety of useful
 *  functions for drawing random numbers from distributions, see
 *  www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html.
 *
 *  Some of the functions are passed a pointer where they should write their
 *  output. The abc_smc function pre-allocates the correct amount of space, so
 *  these functions should not call `malloc` on these pointers. If they
 *  allocate any other memory for local use, it must be freed within the
 *  function.
 *
 *  All the functions can take a user-defined argument, which must be assigned
 *  to the relevant field in the smc_config struct. This is intended to allow
 *  you to adjust the functions based on user input or command line arguments.
 */
typedef struct {
    int nparam; /**< Number of parameters in the model to be fitted */
    int nparticle; /**< Number of particles used to approximate the posterior */
    int nsample; /**< Number of sampled data points per particle */
    int ess_tolerance; /**< ESS below this value triggers resampling */

    double final_epsilon; /**< Tolerance level to end at */
    double final_accept_rate; /**< MCMC acceptance rate to stop at */
    double quality; /**< Between 0 and 1, where 0 is fast and coarse, 1 is slow and accurate */
    double step_tolerance; /**< Tolerance for bisection solution of next epsilon */

    size_t dataset_size; /**< Size of data objects */
    size_t feedback_size; /**< Size of feedback objects */

    /** Perturb a parameter particle.
     *
     *  \param[in] rng       GSL random number generator object, must be used
     *                       for random number generation
     *  \param[in,out] theta Old parameter values, should be overwritten with
     *                       proposed values
     *  \param[in] feedback  Feedback collected from particle population (see
     *                       feedback() ), such as the variance
     *  \param[in] arg       Additional user-defined argument (may be NULL)
     */
    void   (*propose)           (gsl_rng *rng, double *theta, const void *feedback, const void *arg);

    /** Probability density of a proposal.
     *
     *  \param[in] from      Current particle values
     *  \param[in] to        Proposed new particle values
     *  \param[in] feedback  Feedback collected from particle population (see
     *                       feedback())
     *  \param[in] arg       Additional user-defined argument (may be NULL)
     *
     *  \return The probability density of the move from -> to
     */
    double (*proposal_density)  (const double *from, const double *to, const void *feedback, const void *arg);

    /** Simulate a dataset with the given model parameters.
     *
     * \param[in] rng        GSL random number generator, must be used for
     *                       generating any random numbers
     * \param[in] theta      Parameters used to simulate the dataset
     * \param[in] arg        Additional user-defined argument (may be NULL)
     * \param[out] X         The sampled dataset should be stored here
     */
    void   (*sample_dataset)    (gsl_rng *rng, const double *theta, const void *arg, void *X);

    /** Calculate the distance between two datasets.
     *
     * \param[in] a          GSL random number generator, must be used for
     *                       generating any random numbers
     * \param[in] b          Parameters used to simulate the dataset
     * \param[in] arg        Additional user-defined argument (may be NULL)
     *
     * \return The distance between a and b
     */
    double (*distance)          (const void *, const void *, const void *);

    /** Calculate feedback from the population of particles
     *
     * This allows the caller to calculate information about the particles on
     * the fly, in order to modify the proposals. Del Moral et al. (2012) use
     * Gaussian proposals with variance equal to twice the empirical variance
     * of the particles. To achieve this, the feedback function is used to
     * calculate the variance, which is then passed to the propose and
     * proposal_density functions.
     *
     * The population is passed as a single one-dimensional array of length k
     * times n, where k is the number of model parameters and n is the number
     * of particles. The first k elements of the array correspond to the first
     * particle, and the jth parameter of the ith particle is at index k * i +
     * j. 
     *
     * \param[in] theta      Population of parameter particles
     * \param[in] nparticle  Number of particles in the population
     * \param[in] fdbk       Feedback should be stored here
     * \param[in] arg        Additional user-defined argument (may be NULL)
     */
    void   (*feedback)          (const double *theta, int nparticle, void *fdbk, const void *arg);

    /** De-allocate memory for a dataset.
     *
     * \param[in] X          Dataset to free
     */
    void   (*destroy_dataset)   (void *X);

    /** Sample one particle from the prior distribution
     *
     * \param[in] rng        GSL random number generator, must be used for all
     *                       random number generation
     * \param[out] theta     The sampled particle should be copied here.
     * \param[in] arg        Additional user-defined argument (may be NULL)
     */
    void   (*sample_from_prior) (gsl_rng *rng, double *theta, const void *arg);

    /** Prior probability density of a particle
     *
     *  \param[in] theta     Particle values
     *  \param[in] arg       Additional user-defined argument (may be NULL)
     *
     *  \return The prior probability density of theta
     */
    double (*prior_density)     (double *theta, const void *arg);

    void *propose_arg; /**< Extra argument to propose */
    void *proposal_density_arg; /**< Extra argument to proposal_density */
    void *sample_dataset_arg; /**< Extra argument to sample_dataset */
    void *distance_arg; /**< Extra argument to distance */
    void *feedback_arg; /**< Extra argument to feedback */
    void *sample_from_prior_arg; /**< Extra argument to sample_from_prior */
    void *prior_density_arg; /**< Extra argument to prior_density */
} smc_config;

typedef struct {
    int niter;
    double *epsilon;
    double *acceptance_rate;
    double **theta;
    double **W;
} smc_result;

/** Perform ABC-SMC.
 *
 * This implements the adaptive tolerance algorithm from DelMoral et al. 2012.
 * The object returned must be passed to smc_result_free() when you are done
 * with it.
 *
 * \param[in] config control parameters for the algorithm
 * \param[in] seed random seed (if negative, use time)
 * \panam[in] nthread number of threads to use
 * \param[in] data true data
 * \param[in] trace_file file to record particle populations to at each iteration
 * \return an smc_result object with the posterior distribution (theta) and
 * other information about the run
 */
smc_result *abc_smc(const smc_config config, int seed, int nthread, 
                    const void *data, FILE *trace_file);

/** Free an smc_result object.
 *
 * \param[in] r object to free
 */
void smc_result_free(smc_result *r);

#endif
