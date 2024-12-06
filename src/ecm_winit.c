#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "qclib.h"
#include "alloc_utils.h"
#include <string.h>
#include <R_ext/Error.h>

#define BETA 0.05		//Downscaling factor of initial marginal variances
#define MADCONST 1.4826		// normalizing constant for the MAD

static double *__ecm_winit(double *data, int N,
			   int P, int K,
			   int kmed_nrestart,
			   int kmed_itermax, double kmed_tol);

/*
 *  ecm_winit
 *
 *  Wrapper: initializes taus (__ecm_winit) and passes to ecm
 *
 *  Arguments:
 *  - [in] double *data: pointer to data
 *  - [in] int N, P, K: data rows, data cols, num. of clusters
 *  - [in] ecm_*: input passed to ecm function
 *  - [in] kmed_*: input passed to __ecm_winit for kmedian
 *
 *  Returns:
 *  - Return: ECMout (pointer to clustering solution)
 *
 */
ECMout ECM_winit(double *data, int N, int P,
		 int K, double ecm_erc, int ecm_itermax,
		 double ecm_tol, int kmed_nrestart,
		 int kmed_itermax, double kmed_tol)
{

    double *taus = __ecm_winit(data, N, P, K,
			       kmed_nrestart, kmed_itermax,
			       kmed_tol);

    ECMout res;
    if (!taus) {		/* failed init */
	res.info = 100;
	res.iter = 0;
	res.N = N;
	res.P = P;
	res.K = K;
    } else {

	res = ECM(data, N, P, K, taus, ecm_erc, ecm_itermax, ecm_tol);
    }
    return res;
}

/*
 *  __ecm_winit
 *
 *  Initializes the taus weights for ECM algorithm.
 *
 *  Arguments:
 *  - [in] double *data : pointer to data matrix (P x N)
 *  - [in] int N, P     : data dimension (row and cols)
 *  - [in] int K        : number of clusters
 *  - [in] int K        : number of clusters
 *  - [in] int kmed_nrestart: number of times kmedian is run on different initializtion
 *  - [in] double kmed_tol: tolerance for cost function reduction in kemdian
 *  - [in] double kmed_itermax: maximum number of iteration of Lloyd algo in kmedian
 *
 *  Returns:
 *  - [out] double *taus: initial cluster assignment (K x N)
 *
 *  Notes:
 *  - The function calls kmedians (kmedians.c) to initialize the centers, and then
 *    computes diagonal covariances, using scaled down MAD along marginals
 *    Taus are computed as posterior probabilities of points belonging to components
 *    of a mixture of normal distributions with centers and diagonal covariances
 *    computed as described above.
 */
static double *__ecm_winit(double *data, int N,
			   int P, int K,
			   int kmed_nrestart,
			   int kmed_itermax, double kmed_tol)
{

    /* Store log-determinants for checking singular matrices
       (need do be defined before goto statements)
     */
    double logdetS[K];
    memset(logdetS, 0, K * sizeof(double));


    /* initialize table of allocated arrays */
    alloctable *head = NULL;

    double *taus = calloc(K * N, sizeof(double));
    if (taus) {
	alloctable_add(&head, taus, 1, 'A');
    } else {
	warning
	    ("ecm_winit.c (ERR_MALLOC): not able to allocate taus (%d, %d)-vector",
	     K, N);
	goto output;
    }

    double *sumtaus = calloc(N, sizeof(double));
    if (sumtaus) {
	alloctable_add(&head, sumtaus, 0, 'B');
    } else {
	warning
	    ("ecm_winit.c (ERR_MALLOC): not able to allocate sumtaus (%d)-vector",
	     N);
	taus = NULL;
	goto output;
    }

    double *LL = malloc(K * P * sizeof(double));	/* EigenValues */
    if (LL) {
	alloctable_add(&head, LL, 0, 'C');
    } else {
	warning
	    ("ecm_winit.c (ERR_MALLOC): not able to allocate LL (%d, %d)-vector",
	     K, P);
	taus = NULL;
	goto output;
    }




    /* Centers */
    double *MM = kmedians(data, N, P, K, kmed_nrestart, kmed_itermax,
			  kmed_tol);
    if (MM) {
	alloctable_add(&head, MM, 0, 'D');
    } else {
	warning("ecm_winit.c: kmedians failed");
	taus = NULL;
	goto output;
    }



    /* Compute MAD */
    for (int j = 0; j < P; j++) {
	for (int k = 0; k < K; k++) {
	    for (int i = 0; i < N; i++) {
		/* Use first N el.s of taus as a temporary storage */
		*taus++ = fabs(*data++ - *(MM + k * P + j));
	    }
	    data -= N;
	    taus -= N;
	    double lamb_kj =
		BETA * pow(MADCONST * quickselect(taus, N, N / 2), 2.0);
	    *(LL + k * P + j) = lamb_kj;
	    logdetS[k] += log(lamb_kj);
	}
	data += N;
    }
    data -= P * N;

    /* Check singular matrices */
    for (int k = 0; k < K; k++) {
	if (logdetS[k] == -INFINITY) {
	    taus = NULL;
	    goto output;
	}
    }


    /* Compute taus
     *
     * Use the fact that Covariances are diagonal (i.e. Eigenvec = I_p)
     * Use the fact that pi_k = 1/k
     *
     */

    for (int i = 0; i < N; i++)
	*taus++ = 0;
    taus -= N;

    for (int k = 0; k < K; k++) {

	for (int j = 0; j < P; j++) {

	    double mu_kj = *MM++;
	    double lamb_kj = *LL++;

	    /* Squared Mahalanobis Distances */
	    for (int i = 0; i < N; i++)
		*taus++ += pow((*data++ - mu_kj), 2.0) / lamb_kj;
	    taus -= N;
	}
	data -= P * N;

	for (int i = 0; i < N; i++) {
	    *taus = exp(-0.5 * (logdetS[k] + *taus));
	    taus++;
	}

    }
    taus -= K * N;
    MM -= K * P;
    LL -= K * P;

    /* Sumtaus (normalization) */
    for (int k = 0; k < K; k++) {
	for (int i = 0; i < N; i++)
	    *sumtaus++ += *taus++;
	sumtaus -= N;
    }
    taus -= K * N;

    /* Divide taus by sumtaus (normalization) */
    for (int k = 0; k < K; k++) {
	for (int i = 0; i < N; i++)
	    *taus++ /= *sumtaus++;
	sumtaus -= N;
    }
    taus -= K * N;

  output:
    ;
    if (taus)
	alloctable_free_onsuccess(&head);
    else {
	alloctable_free(&head);
    }
    return taus;
}
