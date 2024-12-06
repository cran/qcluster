/* Computes the Smooth Scoring for input: prop, mean, covariances */
#include "qclib.h"
#include "alloc_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <R_ext/Error.h>

/*  score_hard
 *  Hard Scoring with prop, mean and covariances in input.
 *
 *  HS = (nk)^-1 * sum_n max_k(qs(x_i, theta_k))
 *  qs_ik = pi_k - 1/2*log(det(Sigma_k)) - 1/2 [(x_i - mu_k)' inv(Sigma_k) (x_i - mu_k)]
 *
 *  Arguments:
 *  - [in] double *data : pointer to data (a (P,N) matrix)
 *  - [in] double *prop : pointer to vector of proportions (K)
 *  - [in] double *mean : pointer to matrix of centers (K,P)
 *  - [in] double *cov  : pointer to array of covariances (K,P,P)
 *  - [in] int N        : number of observations in data
 *  - [in] int P        : number of dimensions (features) in data
 *  - [in] int K        : number of clusters in solution
 *
 *  Returns double:
 *  the smooth scoring computed on data, using the clustering solution
 *  represented by the triplets: prop, mean, cov
 *
 *  Notes:
 *  qs_ik is computed for all i,k but only the current maximal value is stored
 *  in an array of dimension (N); qs_tmp is used only to facilitate calculations
 *
 *  In computation, we in fact compute 2*QH (to avoid dividing each qs entry by 2)
 *  later, we return the results divided by 2
 */
double score_hard(double *data, int N, int P, int K, double *prop,
		  double *mean, double *cov)
{

    /* Init alloc table */
    alloctable *head = NULL;

    /* Set variables */
    double res = NAN;
    double *qs = malloc(N * sizeof(double));	/* Store quadratic score */
    if (qs)
	alloctable_add(&head, qs, 0, 'A');
    else {
	warning
	    ("score_hard.c (ERR_MALLOC): not able to allocate qs (%d)-vector",
	     N);
	goto output;
    }

    double *Delta = malloc(P * N * sizeof(double));	/* Store mean-centered data */
    if (Delta)
	alloctable_add(&head, Delta, 0, 'B');
    else {
	warning
	    ("score_hard.c (ERR_MALLOC): not able to allocate Delta (%d, %d)-vector",
	     P, N);
	goto output;
    }


    double *covk = malloc(P * P * sizeof(double));	/* Store k-th cov matrix    */
    if (covk)
	alloctable_add(&head, covk, 0, 'C');
    else {
	warning
	    ("score_hard.c (ERR_MALLOC): not able to allocate covk (%d, %d)-vector",
	     P, P);
	goto output;
    }


    /* Compute qs for the k-th triplet */
    for (int k = 0; k < K; k++) {

	/* Delta_k and Cov_k */
	for (int j = 0; j < P; j++) {
	    for (int i = 0; i < N; i++)
		*Delta++ = *data++ - *mean;
	    mean++;

	    for (int i = 0; i <= j; i++)
		*(covk + j * P + i) = *(cov + j * P + i);
	}
	cov += P * P;
	data -= P * N;
	Delta -= P * N;

	/* Compute SMD for k-th cluster */
	double logdet = helper_score_smd(k, N, P, Delta, covk);
	if (isnan(logdet))
	    return NAN;

	/* Compute Sum_k (Z^2) */
	double *qs_tmp = calloc(N, sizeof(double));
	if (!qs_tmp) {
	    warning
		("score_hard.c (ERR_MALLOC): not able to allocate qs_tmp (%d)-vector",
		 N);
	    goto output;
	}

	for (int i = 0; i < P; i++) {
	    for (int j = 0; j < N; j++) {
		*qs_tmp++ += *Delta * *Delta;
		Delta++;
	    }
	    qs_tmp -= N;
	}
	Delta -= P * N;

	/* Fill qs with qs_tmp, if new value are bigger */
	double newval;
	double logprop = 2 * log(*prop++);
	/* Initialize qs, when k==1; else replace old values */
	if (k == 0) {
	    for (int i = 0; i < N; i++)
		*qs++ = logprop - logdet - *qs_tmp++;
	} else {
	    for (int i = 0; i < N; i++) {
		newval = logprop - logdet - *qs_tmp++;
		*qs = newval > *qs ? newval : *qs;
		qs++;
	    }
	}
	qs -= N;
	qs_tmp -= N;
	free(qs_tmp);

    }				/* Close loop over K */

    /* Compute sum(qs) */
    res = 0;
    for (int i = 0; i < N; i++)
	res += *qs++;
    qs -= N;

    res /= (2 * N);

  output:
    ;
    alloctable_free(&head);
    return res;
}
