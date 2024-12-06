/* Computes the Smooth Scoring for input: prop, mean, covariances */
#include "qclib.h"
#include "alloc_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R_ext/Error.h>

/*  score_smooth
 *  Smooth Scoring with prop, mean and covariances in input.
 *
 *  SS = (nk)^-1 * sum_n sum_k tau_ik * qs(x_i, theta_k)
 *  tau_ik = e{qs_ik} / (sum_k e{qs_ik})
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
 *  qs_ik is computed for all i,k and stored in matrix qs of dimension (K,N)
 *  taus_ik analogously --> taus is dimension (K,N)
 *  averages are computed operating on these matrices
 *
 *  In computation, we in fact compute 2*QS (to avoid dividing each qs entry by 2)
 *  later, we return the results divided by 2
 */
double score_smooth(double *data, int N, int P, int K, double *prop,
		    double *mean, double *cov)
{

    /* Init alloc table */
    alloctable *head = NULL;

    /* Set variables */
    double res = NAN;
    double *qs = calloc(K * N, sizeof(double));	/* Store quadratic score */
    if (qs)
	alloctable_add(&head, qs, 0, 'A');
    else {
	warning
	    ("score_smooth.c (ERR_MALLOC): not able to allocate qs (%d)-vector",
	     N);
	goto output;
    }

    double *sumtau = calloc(N, sizeof(double));	/* Store deonminator of tau */
    if (sumtau)
	alloctable_add(&head, sumtau, 0, 'B');
    else {
	warning
	    ("score_smooth.c (ERR_MALLOC): not able to allocate sumtau (%d)-vector",
	     N);
	goto output;
    }

    double *Delta = malloc(P * N * sizeof(double));	/* Store mean-centered data */
    if (Delta)
	alloctable_add(&head, Delta, 0, 'C');
    else {
	warning
	    ("score_smooth.c (ERR_MALLOC): not able to allocate Delta (%d, %d)-vector",
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

	/* Compute Sum_k (Z^2), and store it in k_th row of qs */
	for (int i = 0; i < P; i++) {
	    for (int j = 0; j < N; j++) {
		*qs++ += *Delta * *Delta;
		Delta++;
	    }
	    qs -= N;
	}
	Delta -= P * N;

	/* Compute QS(i,k), and store exp{QS(i,k)}QS(i,k) in qs,
	 * and conteporaneously build sumtau
	 */
	double qs_entry, exp_qs_entry;
	double logprop = 2 * log(*prop);
	for (int i = 0; i < N; i++) {
	    qs_entry = logprop - logdet - *qs;
	    exp_qs_entry = exp(qs_entry);
	    *qs++ = exp_qs_entry * qs_entry;
	    *sumtau++ += exp_qs_entry;
	}
	prop++;
	sumtau -= N;
    }				/* Close loop over K */
    qs -= K * N;

    /* Compute sum(qs / broadcast(sumtau)) */
    res = 0;
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < N; j++) {
	    res += *qs++ / *sumtau++;
	}
	sumtau -= N;
    }
    qs -= K * N;

    res /= (2 * N);

    output:
    ;
    alloctable_free(&head);
    return res;
}
