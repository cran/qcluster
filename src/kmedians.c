#include "qclib.h"
#include "alloc_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R_ext/Error.h>
/*
 *  kmedians
 *
 *  Run kmedians algorithm (i.e. kmeans with L1-distance)
 *
 *  Arguments:
 *  - [in] double *data: pointer to data, transposed format (P X N)
 *  - [in] int N: data number of rows;
 *  - [in] int P: data number of columns;
 *  - [in] int K: number of clusters;
 *  - [in] int random_nrestart: number of times kmedians is run with different random initializ.
 *  - [in] double tol: used to stop the algorithm if cost does not import
 *  - [in] int itermax: maximum number of iterations after which kmedians is stopped
 *
 *  Returns:
 *  - Return: pointer to double (double *MM) of size (K,P); the matrix of optimal centers, with
 *       different clusters' centers stored across rows.
 *
 *  Notes:
 *  1. This function runs a simple version of the kmedian
 *     algorithm, using Lloyd algorithm and random center
 *     initialization.
 *  2. The function is a wrapper to __kmedians, which effectly run the kmedian algorithm once.
 *     The wrapper is used to handle multiple random intializations of the algorithm, and to
 *     return the best one.
 */

static double __kmedians(double *data, int N, int P,
			 int K, int itermax, double tol,
			 double *MM, int *assigned_cluster,
			 double *cluster_marginal, double *l1_dist);

double *kmedians(double *data, int N, int P,
		 int K, int random_nrestart, int itermax, double tol)
{

    /* Init alloc table */
    int success = 1;
    alloctable *head = NULL;

    /*
     * Initialize output
     */
    double *MM = malloc(K * P * sizeof(double));	/* Store cluster centers */
    if (MM)
	alloctable_add(&head, MM, 1, 'A');
    else {
	warning
	    ("kmedians.c (ERR_MALLOC): not able to allocate MM (%d, %d)-array",
	     K, P);
	success = 0;
	goto output;
    }


    /*
     * Init MALLOC-only variables for __kmedians
     */
    double *l1_dist = malloc(N * K * sizeof(double));	/* Store L1-distances from cluster centers */
    if (l1_dist)
	alloctable_add(&head, l1_dist, 0, 'B');
    else {
	warning
	    ("kmedians.c (ERR_MALLOC): not able to allocate l1_dist (%d, %d)-array",
	     N, K);
	success = 0;
	goto output;
    }


    double *cluster_marginal = malloc(N * sizeof(double));
    if (cluster_marginal)
	alloctable_add(&head, cluster_marginal, 0, 'C');
    else {
	warning
	    ("kmedians.c (ERR_MALLOC): not able to allocate cluster_marginal (%d)-array",
	     N);
	success = 0;
	goto output;
    }


    /* int does not fit alloctable */
    int *assigned_cluster = malloc(N * sizeof(int));	/* Store cluster memberships */
    if (!assigned_cluster) {
	warning
	    ("kmedians.c (ERR_MALLOC): not able to allocate assigned_cluster (%d)-array",
	     N);
	success = 0;
	goto output;
    }


    /*
     * Run __kmedians
     */
    double cost[2];
    cost[0] =
	__kmedians(data, N, P, K, itermax, tol, MM, assigned_cluster,
		   cluster_marginal, l1_dist);

    if (random_nrestart > 1) {
	double *MM_new = malloc(K * P * sizeof(double));	/* Store cluster centers */
	if (!MM_new) {
	    warning
		("kmedians.c (ERR_MALLOC): not able to allocate MM_new (%d, %d)-array",
		 K, P);
	    success = 0;
	    goto output;
	}

	/* Re-run kmedians and substitute old one if cost is better */
	for (int i = 1; i < random_nrestart; i++) {
	    cost[1] =
		__kmedians(data, N, P, K, itermax, tol, MM_new,
			   assigned_cluster, cluster_marginal, l1_dist);

	    if (cost[1] < cost[0]) {
		cost[0] = cost[1];
		for (int i = 0; i < K * P; i++)
		    *(MM + i) = *(MM_new + i);
	    }

	}
	free(MM_new);
    }


  output:
    ;
    free(assigned_cluster);
    if (success && (cost[0] < INFINITY)) {
	alloctable_free_onsuccess(&head);
	return MM;
    } else {
	if (!(cost[0] < INFINITY))
	    warning("kmedians.c: not converged");
	alloctable_free(&head);
	return NULL;
    }
}





/*
 *  __kemedians
 *
 *  Run a full kmedian once (initialization, and while loop for Lloyd algo)
 *
 *  Arguments: AS IN KMEDIANS +
 *  - [in,out] double *MM (K,P): the pointer array to store centroids
 *  - [in,out] int *assigned_cluster (N): stores the cluster to which each point is assigned to
 *  - [in,out] double *cluster_marginal (N): auxiliary vector
 *  - [in,out] double *l1_dist (N, K): stores points to clusters distances
 *
 *  Returns: double cost: the overall L1-cost (sum of points' distance to closest center)
 *
 *  Notes:
 *  - All of the additional input this function takes are malloc'ed arrays. Taking them in input,
 *    allows allocating them only once, even if multiple random restarts are required (+efficiency)
 *    Of the additional input, only MM is relevant in output; others get overwritten at each call.
 *
 *  - cluster_marginal is used to store each data feature one at a time, where observations are
 *    ordered by the cluster they are assigned to. This is used in turn to run quickselect on the
 *    relevant portion of the vector (isolating points belonging to each cluster), to compute
 *    marginal medians.
 */
static double __kmedians(double *data, int N,
			 int P, int K, int itermax,
			 double tol, double *MM, int *assigned_cluster,
			 double *cluster_marginal, double *l1_dist)
{


    int cluster_npoints[K];	/* Store cluster number of points */
    for (int k = 0; k < K; k++)
	cluster_npoints[k] = 0;


    double cost[2] = { -1, 0 };	/* Store cost information (old, new) */


    /* Initialization (select K out of N data points) */
    int *pt_ids = Rsample_from(N, K, 0);
    if (!pt_ids)
	return INFINITY;

    for (int k = 0; k < K; k++)
	for (int j = 0; j < P; j++)
	    *(MM + k * P + j) = *(data + *(pt_ids + k) + j * N);
    free(pt_ids);


    int iter = 0;
    while (iter < itermax) {
	iter++;
	/*
	 * STEP-1: FIND CLUSTER ASSIGNMENT FROM NEAREST CENTER
	 */

	/* compute L1 distance (matrix N x K) */
	for (int i = 0; i < N * K; i++)
	    *(l1_dist + i) = 0;
	for (int j = 0; j < P; j++) {
	    for (int i = 0; i < N; i++) {
		for (int k = 0; k < K; k++) {
		    *l1_dist++ += fabs(*data - *(MM + k * P + j));
		}
		data++;
	    }
	    l1_dist -= N * K;
	}
	data -= P * N;

	/* Find nearest center and update cost */
	for (int i = 0; i < N; i++) {
	    double min = *l1_dist++;
	    *assigned_cluster = 0;
	    for (int k = 1; k < K; k++) {
		if (*l1_dist < min) {
		    min = *l1_dist;
		    *assigned_cluster = k;
		}
		l1_dist++;
	    }
	    assigned_cluster++;
	    cost[1] += min;
	}
	l1_dist -= N * K;
	assigned_cluster -= N;



	/*
	 * Check cost and stop while
	 */
	if (fabs(cost[1] - cost[0]) <= tol)
	    break;

	cost[0] = cost[1];
	cost[1] = 0;




	/*
	 * STEP-2: FIND CENTERS FROM CLUSTER ASSIGNMENT
	 * Find marginal median for each cluster, one marginal at a time
	 */

	/* update clusters' info */
	for (int k = 0; k < K; k++)
	    cluster_npoints[k] = 0;
	for (int i = 0; i < N; i++) {
	    cluster_npoints[*(assigned_cluster + i)] += 1;
	}



	/* Find new centers (marginal median of points assigned to each cluster) */
	for (int j = 0; j < P; j++) {

	    int cl_marg_id[K];
	    cl_marg_id[0] = 0;
	    for (int k = 1; k < K; k++)
		cl_marg_id[k] = cluster_npoints[k - 1] + cl_marg_id[k - 1];

	    for (int i = 0; i < N; i++) {
		int k = *(assigned_cluster + i);
		*(cluster_marginal + cl_marg_id[k]) = *data++;
		cl_marg_id[k] += 1;
	    }
	    for (int k = 0; k < K; k++) {
		*(MM + k * P + j) =
		    quickselect((cluster_marginal + cl_marg_id[k] -
				 cluster_npoints[k]),
				cluster_npoints[k],
				cluster_npoints[k] / 2);
	    }

	}

	data -= P * N;

    }

    /* Ensure all clusters have at least 2 points */
    for (int k = 0; k < K; k++)
	if (cluster_npoints[k] < 2)
	    return INFINITY;

    return cost[0];
}
