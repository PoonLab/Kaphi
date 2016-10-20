# include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include "smc.h"
#include "util.h"
#include "stats.h"

#define RESIZE_AMOUNT 100
#define BISECTION_MAX_ITER 10000

/** All the data used for SMC.
 *
 * Because the SMC algorithm is parallelizable, we store all the data in a
 * global structure which is accessed by each thread.
 */
typedef struct {
    const smc_config *config;       /**< configuration parameters and functions */

    const void *data;               /**< input data */
    char *z;                        /**< simulated datasets */
    char *fdbk;                     /**< feedback from current particles */

    double *theta;                  /**< current particles */
    double *new_theta;              /**< storage space for proposed particles */

    double *W;                      /**< current weights */
    double *new_W;                  /**< storage space for proposed new weights */

    double *X;                      /**< distances of simulated datasets */
    double *new_X;                  /**< storage space for proposed new dataset distances */

    double epsilon;                 /**< current tolerance */

    int accept;                     /**< number of accepted proposals */
    int alive;                      /**< number of alive particles */
} smc_workspace;

/** Arguments to the perturb function.
 *
 * One of these objects is created for each thread.
 */
typedef struct {
    int start;        /**< index of first particle (inclusive) */
    int end;          /**< index of last particle (exclusive) */
    int thread_index; /**< thread number */
    gsl_rng *rng;     /**< individual random number generator for this thread */
} thread_data;

/* Use a global workspace instance. */
smc_workspace smc_work;

/* Mutexes for data which will be accessed by many threads. */
pthread_mutex_t smc_accept_mutex;
pthread_mutex_t smc_alive_mutex;

/* Helper functions for SMC. */
void resample(gsl_rng *rng);
double next_epsilon(void);
double epsilon_objfun(double epsilon, void *params);
double ess(const double *W, int n);
void *initialize(void *args);
void *perturb(void *args);

// see Del Moral et al. 2012: An adaptive sequential Monte Carlo method for
// approximate Bayesian computation
smc_result *abc_smc(const smc_config config, int seed, int nthread, 
                    const void *data, FILE *trace_file)
{
    // allocate space for the data in the workspace
    char *z = malloc(config.dataset_size * nthread);
    char *fdbk = malloc(config.feedback_size);
    double *theta = malloc(config.nparticle * config.nparam * sizeof(double));
    double *new_theta = malloc(config.nparticle * config.nparam * sizeof(double));
    double *X = malloc(config.nsample * config.nparticle * sizeof(double));
    double *new_X = malloc(config.nsample * nthread * sizeof(double));
    double *W = malloc(config.nparticle * sizeof(double));
    double *new_W = malloc(config.nparticle * sizeof(double));

    // space for returned values
    double *accept_rate = malloc(RESIZE_AMOUNT * sizeof(double));
    double *epsilons = malloc(RESIZE_AMOUNT * sizeof(double));

    int i, j, niter, status;
    size_t new_size;
    gsl_rng *rng;
    pthread_t *threads = malloc(nthread * sizeof(pthread_t));
    thread_data *thread_args = malloc(nthread * sizeof(thread_data));;
    pthread_attr_t attr;

    smc_result *result = malloc(sizeof(smc_result));
    result->theta = malloc(RESIZE_AMOUNT * sizeof(double*));
    result->W = malloc(RESIZE_AMOUNT * sizeof(double*));

    // initialize pthread things
    pthread_mutex_init(&smc_accept_mutex, NULL);
    pthread_mutex_init(&smc_alive_mutex, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // set up the workspace
    smc_work.config = &config;
    smc_work.data = data;
    smc_work.z = z;
    smc_work.fdbk = fdbk;
    smc_work.theta = theta;
    smc_work.new_theta = new_theta;
    smc_work.W = W;
    smc_work.new_W = new_W;
    smc_work.X = X;
    smc_work.new_X = new_X;
    smc_work.epsilon = DBL_MAX;

    for (i = 0; i < nthread; ++i)
    {
        // we pass start particle (inclusive), end particle (exclusive), and 
        // thread index into each thread
        thread_args[i].start = i * config.nparticle / nthread;
        if (i == nthread - 1) {
            thread_args[i].end = config.nparticle;
        }
        else {
            thread_args[i].end = (i + 1) * config.nparticle / nthread;
        }
        thread_args[i].thread_index = i;

        // also give each of the threads its own random number generator
        thread_args[i].rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(thread_args[i].rng, seed + i + 1);
    }
    rng = set_seed(seed);

    // step 0: sample particles from prior
    smc_work.alive = 0;
    for (i = 0; i < nthread; ++i) {
        status = pthread_create(&threads[i], &attr, initialize, (void *) &thread_args[i]);
    }
    for (i = 0; i < nthread; ++i) {
        status = pthread_join(threads[i], NULL);
    }
    fprintf(stderr, "\n");

    niter = 0;
    printf("iter\tepsilon\tMCMC_accept\n");
    while (smc_work.epsilon != config.final_epsilon)
    {
        printf("%d\t%f\t%f\n", niter, smc_work.epsilon,
               (double) smc_work.accept / (double) smc_work.alive);

        smc_work.accept = 0;
        smc_work.alive = 0;

        // step 1: update epsilon
        smc_work.epsilon = next_epsilon();

        // step 2: resample particles according to their weights
        if (ess(smc_work.W, config.nparticle) < config.ess_tolerance) {
            fprintf(stderr, "ESS = %f, resampling\n", ess(smc_work.W, config.nparticle));
            resample(rng);
        }

        // step 3: perturb particles
        config.feedback(smc_work.theta, config.nparticle, fdbk, config.feedback_arg);
        smc_work.accept = 0;
        smc_work.alive = 0;
        for (i = 0; i < nthread; ++i) {
            status = pthread_create(&threads[i], &attr, perturb, (void *) &thread_args[i]);
        }
        for (i = 0; i < nthread; ++i) {
            status = pthread_join(threads[i], NULL);
        }

        // record everything
        result->theta[niter] = malloc(config.nparticle * config.nparam * sizeof(double));
        memcpy(result->theta[niter], theta, config.nparticle * config.nparam * sizeof(double));
        result->W[niter] = malloc(config.nparticle * sizeof(double));
        memcpy(result->W[niter], W, config.nparticle * sizeof(double));
        epsilons[niter] = smc_work.epsilon;
        accept_rate[niter] = (double) smc_work.accept / (double) smc_work.alive;

        if (trace_file != NULL) {
            for (i = 0; i < config.nparticle; ++i) {
                fprintf(trace_file, "%d\t%f", niter, W[i]);
                for (j = 0; j < config.nparam; ++j) {
                    fprintf(trace_file, "\t%f", theta[i * config.nparam + j]);
                }
                for (j = 0; j < config.nsample; ++j) {
                    fprintf(trace_file, "\t%f", X[i * config.nsample + j]);
                }
                fprintf(trace_file, "\n");
            }
            fflush(trace_file);
        }

        ++niter;

        // allocate more space if necessary
        if (niter % RESIZE_AMOUNT == 0)
        {
            new_size = RESIZE_AMOUNT * (niter / RESIZE_AMOUNT + 1) * sizeof(double);
            epsilons = safe_realloc(epsilons, new_size);
            accept_rate = safe_realloc(accept_rate, new_size);
            new_size = RESIZE_AMOUNT * (niter / RESIZE_AMOUNT + 1) * sizeof(double*);
            result->theta = safe_realloc(result->theta, new_size);
            result->W = safe_realloc(result->W, new_size);
        }

        // if acceptance probability is low enough, we're done
        if (accept_rate[niter-1] <= config.final_accept_rate) {
            smc_work.epsilon = config.final_epsilon;
            break;
        }

    }

    // finally, sample from the estitmated posterior
    resample(rng);
    result->theta[niter] = malloc(config.nparticle * config.nparam * sizeof(double));
    memcpy(result->theta[niter], theta, config.nparticle * config.nparam * sizeof(double));
    result->W[niter] = malloc(config.nparticle * sizeof(double));
    memcpy(result->W[niter], W, config.nparticle * sizeof(double));

    // keep the trace information and final population
    result->niter = niter;
    result->epsilon = epsilons;
    result->acceptance_rate = accept_rate;

    // clean up everything else
    pthread_mutex_destroy(&smc_accept_mutex);
    pthread_mutex_destroy(&smc_alive_mutex);
    pthread_attr_destroy(&attr);
    free(threads);
    free(thread_args);
    for (i = 0; i < nthread; ++i) {
        gsl_rng_free(thread_args[i].rng);
    }
    gsl_rng_free(rng);
    free(z);
    free(fdbk);
    free(theta);
    free(new_theta);
    free(X);
    free(new_X);
    free(W);
    free(new_W);
    return result;
}

void smc_result_free(smc_result *r)
{
    int i;
    free(r->epsilon);
    free(r->acceptance_rate);

    for (i = 0; i <= r->niter; ++i) {
        free(r->theta[i]);
        free(r->W[i]);
    }
    free(r->theta);
    free(r->W);
    free(r);
}

/* Private. */

double next_epsilon(void)
{
    int status, iter = 0;
    double r, x_lo, x_hi;
    double tol = smc_work.config->step_tolerance;
    double prev_epsilon = smc_work.epsilon;
    double alpha = smc_work.config->quality;

    // create a bisection solver
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
    gsl_function F = { .function = &epsilon_objfun, .params = NULL };

    // new epsilon is bounded between 0 and old epsilon
    x_lo = 0;
    x_hi = prev_epsilon;
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    // solve for next epsilon by bisection
    do {
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, tol);
        ++iter;
    }
    while (status == GSL_CONTINUE && iter < BISECTION_MAX_ITER);

    // TODO: handle running out of iterations?
    if (iter == BISECTION_MAX_ITER) {
        fprintf(stderr, "WARNING: hit max iterations solving for next epsilon\n");
        r = smc_work.config->final_epsilon;
    }

    // if epsilon went below the final tolerance or finding the new epsilon
    // failed, then we use the final tolerance
    // we have to call the objective function explicitly to get the new weights
    if (r <= smc_work.config->final_epsilon) {
        r = smc_work.config->final_epsilon;
        epsilon_objfun(r, NULL);
    }

    // use the new weights
    memcpy(smc_work.W, smc_work.new_W, smc_work.config->nparticle * sizeof(double));

    // clean up and return the new epsilon
    gsl_root_fsolver_free(s);
    return r;
}

double epsilon_objfun(double epsilon, void *params)
{
    // reference some of the workspace from global variables for convenience
    double *X = smc_work.X;
    double *W = smc_work.W;
    double *new_W = smc_work.new_W;
    int nparticle = smc_work.config->nparticle;
    int nsample = smc_work.config->nsample;
    double prev_epsilon = smc_work.epsilon;
    double alpha = smc_work.config->quality;

    int i, j; 
    double num, denom, wsum = 0;

    // calculate the new weights (equation 14)
    for (i = 0; i < nparticle; ++i)
    {
        num = 0; 
        denom = 0;
        for (j = 0; j < nsample; ++j)
        {
            num += X[i * nsample + j] < epsilon;
            denom += X[i * nsample + j] < prev_epsilon;
        }
        
        // catch the case when numerator and denominator are both zero
        // TODO: not sure if I should update in that case, or set the weight to zero
        if (num == denom) {
            new_W[i] = W[i];
        }
        else {
            new_W[i] = W[i] * num / denom;
        }
        wsum += new_W[i];
    }

    // normalize the weights so their sum is 1
    for (i = 0; i < nparticle; ++i) {
        new_W[i] /= wsum;
    }

    // if epsilon is too small, the new weights will be undefined
    if (epsilon == 0 || wsum == 0) {
        return -1;
    }

    // equation 12
    // we're trying to minimize ESS(W') - alpha * ESS(W)
    return ess(new_W, nparticle) - alpha * ess(W, nparticle);
}

double ess(const double *W, int n)
{
    int i;
    double sum = 0;
    for (i = 0; i < n; ++i) {
        sum += W[i] * W[i];
    }
    return 1.0 / sum;
}

void resample(gsl_rng *rng)
{
    int i, wcur;
    int nparticle = smc_work.config->nparticle;
    int nparam = smc_work.config->nparam;
    int nsample = smc_work.config->nsample;
    double r, wsum;
    double *W = smc_work.W;
    double *theta = smc_work.theta;
    double *new_theta = smc_work.new_theta;
    double *X = smc_work.X;
    double *new_X = malloc(nparticle * nsample * sizeof(double));

    // sample particles according to their weights
    for (i = 0; i < nparticle; ++i)
    {
        wcur = -1;
        wsum = 0;
        r = gsl_rng_uniform(rng);
        while (r > wsum && wcur < nparticle) {
            wsum += W[++wcur];
        }
        memcpy(&new_theta[i * nparam], &theta[wcur * nparam], nparam * sizeof(double));
        memcpy(&new_X[i * nsample], &X[wcur * nsample], nsample * sizeof(double));
    }

    // reset all the weights
    for (i = 0; i < nparticle; ++i) {
        W[i] = 1.0 / nparticle;
    }

    // use the new particles
    memcpy(smc_work.theta, smc_work.new_theta, nparticle * nparam * sizeof(double));
    memcpy(smc_work.X, new_X, nparticle * nsample * sizeof(double));
    free(new_X);
}

void *initialize(void *args)
{
    int i, j;
    int nparam = smc_work.config->nparam;
    int nparticle = smc_work.config->nparticle;
    int nsample = smc_work.config->nsample;
    thread_data *tdata = (thread_data *) args;
    char *z = &smc_work.z[smc_work.config->dataset_size * tdata->thread_index];
    gsl_rng *rng = tdata->rng;
    int start = tdata->start, end = tdata->end;
    double *particle;

    for (i = start; i < end; ++i)
    {
        particle = &smc_work.theta[i * nparam];
        smc_work.config->sample_from_prior(rng, particle, smc_work.config->sample_from_prior_arg);
        smc_work.W[i] = 1. / nparticle;
        for (j = 0; j < nsample; ++j)
        {
            smc_work.config->sample_dataset(rng, particle, smc_work.config->sample_dataset_arg, z);
            smc_work.X[i * nsample + j] = smc_work.config->distance(z, smc_work.data, smc_work.config->distance_arg);
            smc_work.config->destroy_dataset(z);
        }

        pthread_mutex_lock(&smc_alive_mutex);
        ++smc_work.alive;
        if (smc_work.alive * 10 / nparticle != (smc_work.alive - 1) * 10 / nparticle) {
            fprintf(stderr, "Sampling initial particles (%d%%)\n", smc_work.alive * 100 / nparticle);
        }
        pthread_mutex_unlock(&smc_alive_mutex);
    }
}

void *perturb(void *args)
{
    int i, j;
    double mh_ratio, old_nbhd, new_nbhd;
    double *cur_theta, *prev_theta;
    int nparticle = smc_work.config->nparticle;
    int nparam = smc_work.config->nparam;
    int nsample = smc_work.config->nsample;
    size_t dataset_size = smc_work.config->dataset_size;
    char *fdbk = smc_work.fdbk;
    double epsilon = smc_work.epsilon;

    // get arguments for this thread
    thread_data *tdata = (thread_data *) args;
    int start = tdata->start;
    int end = tdata->end;
    int thread_index = tdata->thread_index;
    gsl_rng *rng = tdata->rng;

    double *W = smc_work.W;
    double *X = smc_work.X;
    double *new_X = &smc_work.new_X[nsample * thread_index];
    char *z = &smc_work.z[dataset_size * thread_index];

    for (i = start; i < end; ++i)
    {
        // ignore dead particles
        if (W[i] == 0)
            continue;

        pthread_mutex_lock(&smc_alive_mutex);
        ++smc_work.alive;
        pthread_mutex_unlock(&smc_alive_mutex);

        cur_theta = &smc_work.new_theta[i * nparam];
        prev_theta = &smc_work.theta[i * nparam];
        memcpy(cur_theta, prev_theta, nparam * sizeof(double));

        // perturb the particle
        smc_work.config->propose(rng, cur_theta, fdbk, smc_work.config->propose_arg);

        // prior ratio
        mh_ratio = smc_work.config->prior_density(cur_theta, smc_work.config->prior_density_arg) /
                   smc_work.config->prior_density(prev_theta, smc_work.config->prior_density_arg);

        if (mh_ratio == 0) {
            continue;
        }

        // proposal ratio
        mh_ratio *= smc_work.config->proposal_density(cur_theta, prev_theta, fdbk, smc_work.config->proposal_density_arg) /
                    smc_work.config->proposal_density(prev_theta, cur_theta, fdbk, smc_work.config->proposal_density_arg);
        if (mh_ratio == 0) {
            continue;
        }

        // sample new datasets
        for (j = 0; j < nsample; ++j)
        {
            smc_work.config->sample_dataset(rng, cur_theta, smc_work.config->sample_dataset_arg, z);
            new_X[j] = smc_work.config->distance(z, smc_work.data, smc_work.config->distance_arg);
            smc_work.config->destroy_dataset(z);
        }

        // SMC approximation to likelihood ratio
        old_nbhd = 0; 
        new_nbhd = 0;
        for (j = 0; j < nsample; ++j) {
            old_nbhd += X[i * nsample + j] < epsilon;
            new_nbhd += new_X[j] < epsilon;
        }
        mh_ratio *= new_nbhd / old_nbhd;

        // accept or reject the proposal
        if (gsl_rng_uniform(rng) < mh_ratio)
        {
            pthread_mutex_lock(&smc_accept_mutex);
            ++smc_work.accept;
            pthread_mutex_unlock(&smc_accept_mutex);

            memcpy(prev_theta, cur_theta, nparam * sizeof(double));
            memcpy(&X[i * nsample], new_X, nsample * sizeof(double));
        }
    }
    return NULL;
}
