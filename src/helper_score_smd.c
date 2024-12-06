#include <stdlib.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Error.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
/*
 *  helper_score_smd
 *
 *  Helps computations of the squared mahalanobis distances efficiently
 *
 *  Arguments:
 *  - [in]     int k        : The cluster being used for SMD
 *  - [in]     int N        : Number of observations
 *  - [in]     int P        : Number of features
 *  - [in,out] double *cov  : Covariance matrix (P,P)
 *  - [in,out] double *Delta: Centered data matrix (P, N)
 *
 *  Returns:
 *  - logdet:  log(det(Sigma)) the log determinant of Cov
 *  - Delta (out): is modified to contain the standardize (non-squared)
 *                 data (P,N): matrix containing (X - mu)/sigma. Summing across P
 *                 gives the Squared-Mahalanobis-Distance vector (N)
 *
 *  NOTE:
 *  - This function is meant to be used by score_smooth and score_hard
 *  - Only lower triangle of cov is referenced
 *  - This function relies on two method Chol2 and Chol4, depending on problem size
 */
static void Chol2(int k, int N, int P,
		  double *Delta, double *cov, double *logdet);
static void Chol4(int k, int N, int P, double *Delta, double *cov);
double helper_score_smd(int k, int N, int P, double *Delta, double *cov)
{

    /* Cholesky Decomposition: Sigma=LL' */
    {
	char uplo = 'U';
	int f77info;
	F77_CALL(dpotrf) (&uplo, &P, cov, &P, &f77info FCONE);
	if (f77info != 0) {
	    warning
		("helper_score_smd (ERR_DPOTRF) on cov_%d exited (info = %d)\nReturning NA\n",
		 k + 1, f77info);
	    return NAN;
	}

    }

    /* Compute Det from Cholesky L.
     * Note: |Sigma| = |L||L'|, computing just |L|
     * we are computing |Sigma|^1/2
     */
    double logdet = 0;
    for (int j = 0; j < P; j++) {
	logdet += 2 * log(*(cov + j * P + j));
    }
    if (fabs(logdet) <= DBL_EPSILON) {
	warning
	    ("helper_score_smd:Chol2 (WARN): %d-th covariance matrix close to singular; returning NA",
	     k + 1);
	return NAN;
    }


    /* Call Chol2 or Chol4 to compute Delta */
    if ((N <= 5000) & (P <= 10))
	Chol2(k, N, P, Delta, cov, &logdet);
    else
	Chol4(k, N, P, Delta, cov);

    return logdet;
}


    /*  Chol2
     *
     *  Computes the SMD as follows
     *  0. Need: Delta %*% Cov^-1 %*% Delta' = rowSums( (Delta) %*% (LL')^-1 * Delta )
     *  1. Use the Cholesky decomposition to compute L
     *  2. Invert the Cholesky factor
     *  3. Compute Z = Delta %*% L^-1
     *
     *  Arguments: Same as calling function + pointer to logdet to handle errors
     *  Returns: Delta (or overwrite logdet on errors)
     *
     */
static void Chol2(int k, int N, int P,
		  double *Delta, double *cov, double *logdet)
{

    /* Compute inverse by inverting Cholesky */
    {
	char uplo = 'U';
	char diag = 'N';
	int f77info;
	F77_CALL(dtrtri) (&uplo, &diag, &P, cov, &P, &f77info FCONE FCONE);

	if (f77info != 0) {
	    warning
		("helper_score_smd:Chol2 (ERR_DPOTRF) on cov_%d exited (info = %d)\n",
		 k + 1, f77info);
	    *logdet = NAN;
	} else {
	    {
		char side = 'R';
		char uplo = 'U';
		char transa = 'N';
		char diag = 'N';
		double alpha = 1;
		F77_CALL(dtrmm) (&side, &uplo, &transa, &diag,
				 &N, &P, &alpha, cov, &P, Delta, &N
				 FCONE FCONE FCONE FCONE);
	    }
	}
    }
}


    /*
     *  Chol4
     *
     *  Computes the SMD as follows
     *  0. Need: Delta %*% Cov^-1 %*% Delta' = (L %*% Delta)(L %*% Delta)' = Z%*%Z'
     *  1. Use Cholesky to decompose Cov: Cov=LL'
     *  2. Use forwardsolve (i.e. solve linear systems with lower-triangular matrix) to solve
     *     L %*% Delta = Z and obtain Z (P, N); Z replaces Delta
     *
     *  Arguments: Same as calling function
     *  Returns: Delta
     *
     */
static void Chol4(int k, int N, int P, double *Delta, double *cov)
{
    /* Solve Linear system Delta = LZ to find SMD */
    {
	char side = 'R';
	char uplo = 'U';
	char transa = 'N';
	char diag = 'N';
	double alpha = 1;
	F77_CALL(dtrsm) (&side, &uplo, &transa, &diag,
			 &N, &P, &alpha, cov, &P, Delta, &N
			 FCONE FCONE FCONE FCONE);
    }
}
