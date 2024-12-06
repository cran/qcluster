#include "qclib.h"
#include "alloc_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <R_ext/BLAS.h>
#include <R_ext/Error.h>
/*
 *  cluster_to_triplets
 *
 *  Computes the (prop,mean,cov) triplets from a (data, cluster) pair
 *
 *
 *  Arguments:
 *  - [in] double *data : pointer to data (a (P,N) matrix)
 *  - [in] int N        : number of observations in data
 *  - [in] int P        : number of dimensions (features) in data
 *  - [in] int K        : number of clusters in solution
 *  - [in] int *cluster : cluster labels
 *
 *  Returns Triplets:
 *  (double) out.prop   : pointer to vector of proportions (K)
 *  (double) out.mean   : pointer to vector of proportions (K, P)
 *  (double) out.cov    : pointer to vector of proportions (K, P ,P)
 *
 */

Triplets cluster_to_triplets(double *data, int N, int P, int K,
			     int *cluster)
{

    int success = 1;
    alloctable *head = NULL;


    double *prop = calloc(K, sizeof(double));
    if (prop)
	alloctable_add(&head, prop, 1, 'A');
    else {
	warning
	    ("cluster_to_triplets.c (ERR_MALLOC): not able to allocate prop (%d)-vector",
	     K);
	success = 0;
	goto output;
    }

    double *mean = calloc(K * P, sizeof(double));
    if (mean)
	alloctable_add(&head, mean, 1, 'B');
    else {
	warning
	    ("cluster_to_triplets.c (ERR_MALLOC): not able to allocate mean (%d, %d)-vector",
	     K, P);
	success = 0;
	goto output;
    }

    double *cov = malloc(K * P * P * sizeof(double));
    if (cov)
	alloctable_add(&head, cov, 1, 'C');
    else {
	warning
	    ("cluster_to_triplets.c (ERR_MALLOC): not able to allocate cov (%d, %d, %d)-vector",
	     K, P, P);
	success = 0;
	goto output;
    }


    /* Start computing number of pts in clusters, and sample means */
    for (int i = 0; i < N; i++) {
	*(prop + *cluster) += 1;
	for (int j = 0; j < P; j++)
	    *(mean + *cluster * P + j) += *(data + i + j * N);
	cluster++;
    }
    cluster -= N;

    for (int k = 0; k < K; k++) {

	for (int j = 0; j < P; j++)
	    *mean++ /= *prop;
	mean -= P;

	int Nk = (int) *prop;
	*prop++ /= N;
	double *points_in_k = malloc(Nk * P * sizeof(double));	/* Subset of points in K; note Col Major (Nk x P) */
	if (!points_in_k) {
	    warning
		("cluster_to_triplets.c (ERR_MALLOC): not able to allocate points_in_k (%d, %d)-vector",
		 Nk, P);
	    success = 0;
	    goto output;
	}

	for (int i = 0; i < N; i++) {
	    if (*cluster == k) {
		for (int j = 0; j < P; j++)
		    *points_in_k++ = *(data + i + j * N) - *(mean + j);	/* centering the data */
	    }
	    cluster++;
	}
	mean += P;
	cluster -= N;
	points_in_k -= Nk * P;

	/* Compute sample covariance */
	{
	    char uplo = 'U';
	    char trans = 'N';
	    double alpha = 1.0 / (Nk - 1);
	    double beta = 0;
	    F77_CALL(dsyrk) (&uplo, &trans, &P, &Nk,
			     &alpha, points_in_k, &P, &beta, cov, &P
			     FCONE FCONE);
	}

	for (int i = 0; i < P; i++)
	    for (int j = 0; j < i; j++)
		*(cov + i + j * P) = *(cov + i * P + j);
	cov += P * P;
	free(points_in_k);
    }				/* Close loop over K */


    prop -= K;
    mean -= K * P;
    cov -= K * P * P;

  output:
    ;
    Triplets out;
    if (success) {
	out.prop = prop;
	out.mean = mean;
	out.cov = cov;
	alloctable_free_onsuccess(&head);
    } else {
	out.prop = NULL;
	out.mean = NULL;
	out.cov = NULL;
	alloctable_free(&head);
    }

    return out;
}
