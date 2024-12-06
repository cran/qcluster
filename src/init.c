#include "qclib.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include <Rinternals.h>

/* Wrapper for ECM (for help in input argument see calling functions from R) */
SEXP ECM_C(SEXP data, SEXP nrow, SEXP ncol, SEXP ncluster, SEXP kmed_init,
	   SEXP ecm_init, SEXP ecm_erc, SEXP ecm_maxiter, SEXP ecm_tol,
	   SEXP kmed_nstart, SEXP kmed_maxiter,
	   SEXP kmed_tol, SEXP save_cluster, SEXP save_params,
	   SEXP save_taus)
{

    /* Perparing variables from R */
    int N = asInteger(nrow);
    int P = asInteger(ncol);
    int K = asInteger(ncluster);
    int init = asInteger(kmed_init);

    /* Obtain pointer to data location */
    double *dt = REAL(data);

    ECMout ans;
    if (init == 1) {
	/* Ignore ecm_init and use kmedian for initialization */
	ans = ECM_winit(dt, N, P, K,
			asReal(ecm_erc), asInteger(ecm_maxiter),
			asReal(ecm_tol), asInteger(kmed_nstart),
			asInteger(kmed_maxiter), asReal(kmed_tol));
    } else {
	/* Use ecm_init as initial weights (need to reallocate weight matrix
	 * R will otherwise complain about double freeing the taus vector
	 */
	double *ini = REAL(ecm_init);
	double *taus_init = malloc(K * N * sizeof(double));
	if (!taus_init) {
	    warning
		("init.c (ERR_MALLOC): not able to allocate double taus_init (%d, %d)-array",
		 K, N);
	    ans.info = 0;
	    ans.iter = 0;
	    ans.N = N;
	    ans.P = P;
	    ans.K = K;
	    goto output;
	}

	for (int i = 0; i < K * N; i++)
	    *taus_init++ = *ini++;
	taus_init -= K * N;
	ans =
	    ECM(dt, N, P, K, taus_init, asReal(ecm_erc),
		asInteger(ecm_maxiter), asReal(ecm_tol));
    }

  output:
    ;
    /* Prepare output (a SEXP list) */
    int prot = 0, nout;
    SEXP output;
    nout = ans.info >= 200 ? 10 : 5;
    output = PROTECT(allocVector(VECSXP, nout));
    prot++;

    /* info */
    SEXP info = PROTECT(allocVector(INTSXP, 1));
    prot++;
    INTEGER(info)[0] = ans.info;
    SET_VECTOR_ELT(output, 0, info);

    /* iter */
    SEXP iter = PROTECT(allocVector(INTSXP, 1));
    prot++;
    INTEGER(iter)[0] = ans.iter;
    SET_VECTOR_ELT(output, 1, iter);

    /* N, P, K */
    SEXP rN = PROTECT(allocVector(INTSXP, 1));
    prot++;
    INTEGER(rN)[0] = ans.N;
    SET_VECTOR_ELT(output, 2, rN);

    SEXP rP = PROTECT(allocVector(INTSXP, 1));
    prot++;
    INTEGER(rP)[0] = ans.P;
    SET_VECTOR_ELT(output, 3, rP);

    SEXP rK = PROTECT(allocVector(INTSXP, 1));
    prot++;
    INTEGER(rK)[0] = ans.K;
    SET_VECTOR_ELT(output, 4, rK);

    if (ans.info >= 200) {

	/* eloglik */
	SEXP eloglik = PROTECT(allocVector(REALSXP, 1));
	prot++;
	REAL(eloglik)[0] = ans.eloglik;
	SET_VECTOR_ELT(output, 5, eloglik);

	/* Cluster sizes */
	SEXP size = PROTECT(allocVector(INTSXP, K));
	prot++;
	int *psize = INTEGER(size);
	for (int i = 0; i < K; i++)
	    *psize++ = *ans.size++;
	ans.size -= K;
	SET_VECTOR_ELT(output, 6, size);
	free(ans.size);

	/* cluster */
	if (INTEGER(save_cluster)[0]) {
	    SEXP cluster = PROTECT(allocVector(INTSXP, N));
	    prot++;
	    int *pcluster = INTEGER(cluster);
	    for (int i = 0; i < N; i++)
		*pcluster++ = *ans.cluster++;
	    ans.cluster -= N;
	    SET_VECTOR_ELT(output, 7, cluster);
	}
	free(ans.cluster);

	/* Taus */
	if (INTEGER(save_taus)[0]) {
	    SEXP taus = PROTECT(allocMatrix(REALSXP, N, K));
	    prot++;
	    double *ptaus = REAL(taus);
	    for (int i = 0; i < N * K; i++)
		*ptaus++ = *ans.tau++;
	    ans.tau -= K * N;
	    SET_VECTOR_ELT(output, 8, taus);
	}
	free(ans.tau);

	/* Params */
	if (INTEGER(save_params)[0]) {
	    SEXP params = PROTECT(allocVector(VECSXP, 3));
	    prot++;
	    {
		SEXP prop = PROTECT(allocVector(REALSXP, K));
		prot++;
		SEXP mean = PROTECT(allocMatrix(REALSXP, P, K));
		prot++;
		SEXP cov = PROTECT(allocVector(REALSXP, K * P * P));	/* Reshaped in R */
		prot++;

		double *pparam;

		pparam = REAL(prop);
		for (int i = 0; i < K; i++)
		    *pparam++ = *ans.params.prop++;
		ans.params.prop -= K;
		free(ans.params.prop);
		SET_VECTOR_ELT(params, 0, prop);

		pparam = REAL(mean);
		for (int i = 0; i < K * P; i++)
		    *pparam++ = *ans.params.mean++;
		ans.params.mean -= K * P;
		free(ans.params.mean);
		SET_VECTOR_ELT(params, 1, mean);

		pparam = REAL(cov);
		for (int i = 0; i < P * P * K; i++)
		    *pparam++ = *ans.params.cov++;
		ans.params.cov -= K * P * P;
		free(ans.params.cov);
		SET_VECTOR_ELT(params, 2, cov);
	    }
	    SET_VECTOR_ELT(output, 9, params);
	} else {
	    free(ans.params.prop);
	    free(ans.params.mean);
	    free(ans.params.cov);
	}

    }

    /* NAMES */
    UNPROTECT(prot);

    return output;
}


/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < 3; i++) {
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    }
    return elmt;
}


/* Scoring function; wrapper for score_*.c (for help in input argument see calling functions from R) */
SEXP SCORE_C(SEXP type, SEXP data, SEXP nrow, SEXP ncol, SEXP ncluster,
	     SEXP params)
{

    /* prepare input */
    int N = asInteger(nrow);
    int P = asInteger(ncol);
    int K = asInteger(ncluster);
    int in_score_type = asInteger(type);
    double *dt = REAL(data);

    int prot = 0;
    /* prepare output */
    SEXP output;
    output = PROTECT(allocVector(REALSXP, 2));
    prot++;

    double hard, smooth;



    /* Calling score functions */
    double *prop = REAL(getListElement(params, "proportion"));
    double *mean = REAL(getListElement(params, "mean"));
    double *cov = REAL(getListElement(params, "cov"));

    if (in_score_type == 2) {
	hard = score_hard(dt, N, P, K, prop, mean, cov);
	smooth = score_smooth(dt, N, P, K, prop, mean, cov);
    } else if (in_score_type == 1) {
	hard = R_NaN;
	smooth = score_smooth(dt, N, P, K, prop, mean, cov);
    } else {
	hard = score_hard(dt, N, P, K, prop, mean, cov);
	smooth = R_NaN;
    }



    REAL(output)[0] = hard;
    REAL(output)[1] = smooth;


    UNPROTECT(prot);
    return (output);
}

/* Reutrn clusters' params triplets from data and labels */
SEXP TRIPLETS_C(SEXP data, SEXP nrow, SEXP ncol, SEXP ncluster,
		SEXP cluster)
{

    /* Handle input */
    int N = asInteger(nrow);
    int P = asInteger(ncol);
    int K = asInteger(ncluster);
    double *dt = REAL(data);
    int *labels = INTEGER(cluster);

    Triplets params = cluster_to_triplets(dt, N, P, K, labels);

    /* Handle output (SEXP list) */
    SEXP output;
    if (!params.prop || !params.mean || !params.cov) {
	output = R_NilValue;
    } else {
	int prot = 0;
	output = PROTECT(allocVector(VECSXP, 3));
	prot++;

	SEXP prop = PROTECT(allocVector(REALSXP, K));
	prot++;
	SEXP mean = PROTECT(allocMatrix(REALSXP, P, K));
	prot++;
	SEXP cov = PROTECT(allocVector(REALSXP, K * P * P));	/* Reshaped in R */
	prot++;

	/* Copy to output */
	double *pparam;

	pparam = REAL(prop);
	for (int i = 0; i < K; i++)
	    *pparam++ = *params.prop++;
	params.prop -= K;

	pparam = REAL(mean);
	for (int i = 0; i < K * P; i++)
	    *pparam++ = *params.mean++;
	params.mean -= K * P;

	pparam = REAL(cov);
	for (int i = 0; i < P * P * K; i++)
	    *pparam++ = *params.cov++;
	params.cov -= K * P * P;

	SET_VECTOR_ELT(output, 0, prop);
	SET_VECTOR_ELT(output, 1, mean);
	SET_VECTOR_ELT(output, 2, cov);

	free_Triplets(params);
	UNPROTECT(prot);
    }
    return (output);
}

    /* Registering methods */
static const R_CallMethodDef callMethods[] = {
    { "ECM_C", (DL_FUNC) & ECM_C, 15 },
    { "SCORE_C", (DL_FUNC) & SCORE_C, 6 },
    { "TRIPLETS_C", (DL_FUNC) & TRIPLETS_C, 5 },
    { NULL, NULL, 0 }
};

void R_init_qcluster(DllInfo * dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
