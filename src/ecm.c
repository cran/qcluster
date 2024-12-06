#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
#include "qclib.h"
#include "alloc_utils.h"
#include <string.h>

#define PI 3.14159265358979323846

/* ECM
 * Runs the ECM algorithm on data; returns prop, mean, covs
 * and clustering
 *
 * Arguments:
 * - [in]     data    : pointer to data matrix (P x N)
 * - [in]     N, P    : data dimension (row and cols)
 * - [in,out] taus    : initial cluster assignment (K x N)
 * - [in]     K       : number of clusters
 * - [in]     erc     : the value of eigenratio constraint
 * - [in]     itermax : maximum iterations for EM
 * - [in]     tol     : EM tolerace to assess convergence
 *
 *
 * Output:
 * - ECMout: struct containing several ECM objects
 *
 *
 * Notes:
 * flags (embedded in "info")
 *     0 = noflag
 *     1 = numerically degenerate tau cannot be prevented
 *     2 = successfully enforced the ERC at least once
 *     3 = 1 + 2
 *
 * info: 3 digits: 1st digit (hundreds) is for code
 *   -200: used to signal DSYEV fail (returns NA)
 *   -100: used to signal MALLOC fail (returns NA)
 *    10X: no better than initial
 *    20X: iter.max reached
 *    30X: convergence within iter.max
 *     (X): where X is from flags.
 *
 *    e.g.: 101 = code: No better than initial; flag: numerically degenerate
 *
 * Code objects dimensions:
 *    taus, taus_new: (K x N)
 *    MM, LL        : (K x P)
 *    data, Delta, Xw: (P x N)
 *    VV: (K x P x P)
 *    sumtau, sumtau_new: (K)
 */
ECMout ECM(double *data, int N, int P,
	   int K, double *taus, double erc, int itermax, double tol)
{

    /* Table of allocated vectors */
    alloctable *head = NULL;
    alloctable_add(&head, taus, 1, '9');

    /* Define variables */
    int info = 0;
    double *sumtau = calloc(K, sizeof(double));	/* sumtau (alloc) */
    alloctable *first_head = NULL;
    if (sumtau) {
	alloctable_add(&head, sumtau, 1, 'A');
	first_head = head;
    } else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate sumtau (%d)-vector",
	     K);
	info = -100;
	goto output;
    }
    const double GAUSSCOST = pow(2.0 * PI, -(double) P / 2.0);

    /* Initialize sumtau */
    for (int i = 1; i <= K * N; i++) {
	*sumtau += *taus++;
	if ((i % N) == 0) {
	    sumtau++;
	}
    }
    taus -= N * K;		/* Reset pointer */
    sumtau -= K;

    int dowhile = 1;		/* replaces "stop" variable in R */
    int iter = 0;

    /* Cluster means (don't need to be initialized at 0) */
    double *MM = malloc(K * P * sizeof(double));	/* MM (alloc) */
    if (MM)
	alloctable_add(&head, MM, 1, 'B');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate MM (%d,%d)-vector",
	     K, P);
	info = -100;
	goto output;
    }

    double *cov = NULL;

    /* Prepare Delta and Xpp (reused over clusters) */
    double *Delta = malloc(P * N * sizeof(double));
    if (Delta)
	alloctable_add(&head, Delta, 0, 'C');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate Delta (%d, %d)-vector",
	     P, N);
	info = -100;
	goto output;
    }

    double *SS = malloc(P * P * sizeof(double));
    if (SS)
	alloctable_add(&head, SS, 0, 'D');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate SS (%d, %d)-vector",
	     P, P);
	info = -100;
	goto output;
    }

    /* These variables need to be initialized at 0 in each iteration of while loop (allocate once) */
    double *VV = malloc(K * P * P * sizeof(double));	/* VV (alloc) */
    if (VV)
	alloctable_add(&head, VV, 1, 'E');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate VV (%d, %d, %d)-vector",
	     K, P, P);
	info = -100;
	goto output;
    }

    double *LL = malloc(K * P * sizeof(double));	/* LL (alloc) */
    if (LL)
	alloctable_add(&head, LL, 0, 'F');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate LL (%d, %d)-vector",
	     K, P);
	info = -100;
	goto output;
    }

    double *psi = malloc(N * sizeof(double));	/* psi (alloc) */
    if (psi)
	alloctable_add(&head, psi, 0, 'G');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate psi (%d)-vector",
	     N);
	info = -100;
	goto output;
    }

    /* Temporary variables */
    double *Xw = malloc(P * N * sizeof(double));	/* Xw (alloc) */
    if (Xw)
	alloctable_add(&head, Xw, 0, 'H');
    else {
	warning
	    ("ecm.c (ERR_MALLOC): not able to allocate Xw (%d, %d)-vector",
	     P, N);
	info = -100;
	goto output;
    }

    /* lwork for DSYEV (lapack) */
    int lwork = -1;
    int f77info = 0;

    /* EigenRatio constraint variables */
    double Lmin, Lmax, EigenRatio;

    /* Log-likelihood  */
    double eloglik = -INFINITY;

    while (dowhile) {

	/* These variables need to be initialized to 0 at each iteration */
	memset(VV, 0, K * P * P * sizeof(double));
	memset(LL, 0, K * P * sizeof(double));
	memset(psi, 0, N * sizeof(double));
	double *taus_new = calloc(K * N, sizeof(double));	/* taus_new (alloc) */
	if (taus_new)
	    alloctable_add(&head, taus_new, 0, 'I');
	else {
	    warning
		("ecm.c (ERR_MALLOC): not able to allocate taus_new (%d, %d)-vector",
		 K, N);
	    info = -100;
	    goto output;
	}

	double *sumtau_new = calloc(K, sizeof(double));	/* sumtau_new (alloc) */
	if (sumtau_new)
	    alloctable_add(&head, sumtau_new, 0, 'J');
	else {
	    warning
		("ecm.c (ERR_MALLOC): not able to allocate sumtau_new (%d)-vector",
		 K);
	    info = -100;
	    goto output;
	}


	for (int k = 0; k < K; k++) {

	    /*  Prepare K-th row of MM (and initialize LL and VV) */
	    for (int j = 0; j < P; j++) {

		/*  Prepare K-th row of MM */
		double Mentry = 0;	/* Mentry */
		for (int i = 0; i < N; i++) {
		    Mentry += *data++ * *taus++;
		}
		taus -= N;
		*MM++ = Mentry / *sumtau;
	    }
	    data -= P * N;
	    MM -= P;

	    for (int j = 0; j < P; j++) {
		for (int i = 0; i < N; i++) {
		    *Delta = *data++ - *MM;
		    *Xw++ = *Delta++ * sqrt(*taus++ / *sumtau);
		}
		MM++;
		taus -= N;
	    }
	    Xw -= P * N;

	    /* Compute  crossprod */
	    {
		char uplo = 'U';
		char trans = 'T';
		double alpha = 1;
		double beta = 1;
		F77_CALL(dsyrk) (&uplo, &trans, &P, &N,
				 &alpha, Xw, &N, &beta, VV, &P
				 FCONE FCONE);
	    }

	    /* Compute Eigenval Eigenvec */
	    if (f77info == 0) {
		char jobz = 'V';
		char uplo = 'U';
		double size_work;
		/* first call to query optimal lwork */
		if (lwork == -1) {
		    F77_CALL(dsyev) (&jobz, &uplo, &P, VV, &P, LL,
				     &size_work, &lwork, &f77info
				     FCONE FCONE);
		    lwork = (int) size_work;
		}

		double *work = calloc(lwork, sizeof(double));
		if (!work) {
		    warning
			("ecm.c (ERR_MALLOC): not able to allocate work (%d)-vector",
			 lwork);
		    info = -100;
		    goto output;
		}

		F77_CALL(dsyev) (&jobz, &uplo, &P, VV, &P, LL, work,
				 &lwork, &f77info FCONE FCONE);
		free(work);
		if (f77info != 0) {
		    warning("ecm.c (ERR_DSYEV): info = %d\n", f77info);
		    info = -200;
		    goto output;
		}
	    }

	    /* Set pointer for next iteration */
	    taus += N;		/* Move to next block */
	    sumtau += 1;	/* Move to next block */
	    data -= P * N;	/* Return to beginning of data */
	    VV += P * P;	/* Move to next block */
	    LL += P;		/* Move to next block */
	    Delta -= P * N;	/* overwrite */

	}			/*  Close first loop over K */
	/* Reset Pointers */
	LL -= K * P;
	VV -= K * P * P;
	taus -= K * N;
	sumtau -= K;
	MM -= K * P;

	/* Implementing EigenRatio constraint */
	if (isfinite(erc)) {
	    Lmin = *LL;
	    Lmax = *(LL + P - 1);
	    for (int k = 1; k < K; k++) {
		if (Lmin > *(LL + P * k))
		    Lmin = *(LL + P * k);
		if (Lmax < *(LL + P - 1 + P * k))
		    Lmax = *(LL + P - 1 + P * k);
	    }
	    EigenRatio = Lmax / Lmin;
	    if ((EigenRatio > erc) || !isfinite(EigenRatio) || Lmin <= 0
		|| !isfinite(Lmax)) {
		/*  if flag contains 2 do not add it */
		{
		    int flag = info % 100;
		    info += flag > 1 ? 0 : 2;
		}
		GssERC_v2(LL, erc, Lmin, Lmax, sumtau, P, K);
	    }
	}


	for (int k = 0; k < K; k++) {
	    /* Compute Sigma inverse using the formula
	     * SSinv = EigVec %*% EigVal^-1 %*% EigVec^T
	     * We use the decomposition:
	     * SSinv = AA', where A = EigVec %*% EigVal^-1/2
	     * After the following loop, SS = t(VV[,,j] %*% diag(1/sqrt(L[,j]))
	     */
	    for (int i = 0; i < P * P; i++)
		*SS++ = *VV++ / pow(*(LL + i / P), 0.5);
	    SS -= P * P;
	    /*  LL += P; This prepares for next iter, but is auto-handled in next for loop (phi) */

	    /* NOTE: blas accumulates values in SSinv;
	     * need to be redefined each time */
	    double *SSinv = calloc(P * P, sizeof(double));	/* SSinv (alloc) */
	    if (SSinv)
		alloctable_add(&head, SSinv, 0, 'L');
	    else {
		warning
		    ("ecm.c (ERR_MALLOC): not able to allocate SSinv (%d, %d)-vector",
		     P, P);
		info = -100;
		goto output;
	    }

	    {
		char uplo = 'L';
		char trans = 'N';
		double alpha = 1;
		double beta = 1;
		F77_CALL(dsyrk) (&uplo, &trans, &P, &P,
				 &alpha, SS, &P, &beta, SSinv,
				 &P FCONE FCONE);
	    }


	    /* Squared Mahalanobis distances */
	    /*  Delta (create column wise) */
	    for (int j = 0; j < P; j++) {
		for (int i = 0; i < N; i++)
		    *Delta++ = *data++ - *MM;
		MM++;
	    }
	    Delta -= P * N;
	    data -= P * N;

	    /*  Delta %*% SSinv */
	    double *DS = malloc(P * N * sizeof(double));
	    if (DS)
		alloctable_add(&head, DS, 0, 'M');
	    else {
		warning
		    ("ecm.c (ERR_MALLOC): not able to allocate DS (%d, %d)-vector",
		     P, N);
		info = -100;
		goto output;
	    }

	    /* NOTE: in Fortran Delta is seen as NxP, thus we write: Delta%*%SSinv
	       the matrix DS will be created in colum-mode in fortran (NxP),
	       so it must be treated in C as PxN (as well as Delta)
	     */
	    {
		char side = 'R';
		char uplo = 'L';
		double alpha = 1;
		double beta = 0;
		F77_CALL(dsymm) (&side, &uplo, &N, &P, &alpha, SSinv, &P,
				 Delta, &N, &beta, DS, &N FCONE FCONE);
	    }

	    /*  (Delta %*% SSinv) * Delta */
	    for (int j = 0; j < P; j++) {
		for (int i = 0; i < N; i++) {
		    *Delta = *Delta * *DS++;
		    Delta++;
		}
	    }
	    Delta -= P * N;
	    DS -= P * N;
	    alloctable_free_last(&head);
	    alloctable_free_last(&head);

	    /*  Cluster multiplier */
	    double cluster_multiplier = 1;
	    for (int j = 0; j < P; j++)
		cluster_multiplier *= *LL++;
	    cluster_multiplier =
		*sumtau++ / N * GAUSSCOST * pow(cluster_multiplier, -0.5);

	    /*  Compute Phi's (store them in taus_new) */
	    for (int i = 0; i < N; i++) {
		for (int j = 0; j < P; j++)
		    *taus_new += *(Delta + N * j + i);
		*taus_new = cluster_multiplier * exp(-0.5 * *taus_new);
		taus_new++;
	    }
	}			/*  Close second loop over clusters K */

	/* Reset pointers */
	LL -= K * P;
	VV -= K * P * P;
	MM -= K * P;
	taus_new -= K * N;
	sumtau -= K;

	/* Compute new weights */
	/* psi */
	for (int k = 0; k < K; k++) {
	    for (int i = 0; i < N; i++)
		*psi++ += *taus_new++;
	    psi -= N;
	}
	taus_new -= K * N;

	/* taus_new (with posterior check sub-step) */
	for (int k = 0; k < K; k++) {
	    for (int i = 0; i < N; i++) {
		*taus_new /= *psi++;
		*sumtau_new += *taus_new;
		if (!isfinite(*taus_new)) {
		    /*  add flag 1 and set info to 100 (=em.failed) */
		    info += 100 + info % 100 + 1;
		    goto output;
		}
		taus_new++;
	    }
	    if (*sumtau_new == 0) {
		/*  add flag 1 and set info to 100 (=em.failed) */
		info = 100 + info % 100 + 1;
		goto output;
	    }
	    sumtau_new++;
	    psi -= N;
	}
	sumtau_new -= K;
	taus_new -= K * N;

	/* Posterior check: substeps 2a and 2b in JASA */
	double *tmp = calloc(N, sizeof(double));
	if (!tmp) {
	    warning
		("ecm.c (ERR_MALLOC): not able to allocate tmp (%d)-vector",
		 N);
	    info = -100;
	    goto output;
	}

	for (int k = 0; k < K; k++) {
	    for (int i = 0; i < N; i++) {
		*tmp++ += *taus_new++;
	    }
	    tmp -= N;
	}
	taus_new -= K * N;

	/*  Find zeros of tmp; use first positions of tmp to store their locations */
	int n_tmp = 0;
	double *i_tmp = tmp;
	for (int i = 0; i < N; i++) {
	    if (*tmp++ == 0) {
		*i_tmp++ = i;
		n_tmp++;
	    }
	}
	i_tmp -= n_tmp;

	if (n_tmp > 0) {

	    int cl_tmp = 0;
	    double val1 = INFINITY, val2 = 0;
	    for (int i = 0; i < n_tmp; i++) {
		for (int k = 0; k < K; k++) {
		    *(taus_new + (int) *(i_tmp + i) + N * k) = 0;

		    for (int j = 0; j < P; j++) {
			double x_tmp =
			    *(data + (int) *(i_tmp + i) + N * j) -
			    *(MM + j + k * P);
			val2 += x_tmp * x_tmp;
		    }

		    if (val2 < val1)
			cl_tmp = k;
		    val1 = val2;
		    val2 = 0;
		}

		for (int k = 0; k < K; k++)
		    *(taus_new + (int) *(i_tmp + i) + cl_tmp * N) = 1;
	    }
	    for (int k = 0; k < K; k++) {
		*sumtau_new = 0;
		for (int i = 0; i < N; i++)
		    *sumtau_new += *taus_new++;
		sumtau_new++;
	    }
	    sumtau_new -= K;
	    taus_new -= K * N;
	}
	tmp -= N;
	free(tmp);


	iter++;
	/*  loglik check */
	double eloglik_new = 0, eloglik_delta;
	for (int i = 0; i < N; i++) {
	    eloglik_new += log(*psi++);
	}
	psi -= N;

	eloglik_delta = fabs(eloglik_new - eloglik);
	if ((iter >= itermax)
	    || ((eloglik_delta) && (eloglik_delta < tol))) {
	    dowhile = 0;

	    /*  Prepare Covariance in output */
	    for (int k = 0; k < K; k++) {
		for (int i = 0; i < P * P; i++) {
		    *SS++ = *VV * pow(*(LL + i / P), 0.5);
		    *VV++ = 0;
		}
		SS -= P * P;
		VV -= P * P;
		{
		    char uplo = 'L';
		    char trans = 'N';
		    double alpha = 1;
		    double beta = 1;
		    F77_CALL(dsyrk) (&uplo, &trans, &P, &P,
				     &alpha, SS, &P, &beta, VV, &P
				     FCONE FCONE);
		}

		for (int i = 0; i < P; i++) {
		    for (int j = 0; j < P; j++) {
			*(VV + i + j * P) = *(VV + i * P + j);
		    }
		}

		VV += P * P;
		LL += P;
	    }
	    VV -= K * P * P;
	    LL -= K * P;
	    cov = VV;

	}

	eloglik = eloglik_new;


	alloctable *curr;

	/* Sumtau and taus are treated differently: they exchange dynamically
	   with *_new versions. Thus, the *_new versions are stored in the alloctable
	   only to avoid handle alloc failures, but are removed and handled externally
	 */
	free(sumtau);
	if (first_head) {	/* avoid double free in case of errors */
	    first_head->allocated = NULL;
	    first_head = first_head->next;
	}

	sumtau = sumtau_new;
	curr = head;
	head = head->next;
	free(curr);

	free(taus);
	if (first_head) {	/* avoid double free in case of errors */
	    first_head->allocated = NULL;
	    first_head = NULL;
	}

	taus = taus_new;
	curr = head;
	head = head->next;
	free(curr);
    }				/*  close while loop */


    /* Assign solution code (Format flags not needed) */
    if (iter == itermax) {
	info = 200 + info % 100;
    } else if (iter < itermax) {
	info = 300 + info % 100;
    }

  output:
    ;
    /* Output */
    ECMout res;
    res.info = info;
    res.iter = iter;
    res.N = N;
    res.P = P;
    res.K = K;
    if (info < 200) {
	alloctable_free(&head);
	res.cluster = NULL;
	res.size = NULL;
	res.tau = NULL;
	res.params.mean = NULL;
	res.params.cov = NULL;
	res.params.prop = NULL;
    } else {
	alloctable_free_onsuccess(&head);
	res.eloglik = eloglik / N;
	res.params.mean = MM;
	res.params.cov = cov;
	res.tau = taus;

	res.params.prop = sumtau;
	for (int k = 0; k < K; k++) {
	    *sumtau++ /= N;
	}

	/* At this point, should have enough memory to allocate K + N doubles */
	res.cluster = calloc(N, sizeof(int));
	res.size = calloc(K, sizeof(int));

	for (int i = 0; i < N; i++) {
	    for (int k = 0; k < K; k++) {
		*res.cluster =
		    *(taus + i + *res.cluster * N) >
		    *(taus + i + k * N) ? *res.cluster : k;
	    }
	    *(res.size + *res.cluster) += 1;
	    res.cluster++;
	}
	res.cluster -= N;
    }

    return res;
}
