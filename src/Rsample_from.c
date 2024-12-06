#include <stdlib.h>
#include <stdio.h>
#include <R_ext/Random.h>
#include <R_ext/Error.h>
/*
 *  Rsample_from
 *
 *  Description: return K integers from a range [0,N)
 *
 *  Arguments:
 *  - [in] int K: integers to estract
 *  - [in] int N: upper bound (0 to N, in step 1)
 *  - [in] int replace: boolean (0 sample without replacement, 1 sample with replacement)
 *
 *  Returns:
 *  - Return: integer *ret: the K sampled numbers (dim: K)
 *
 *  Notes: this function uses R libraries. Random seed in handled by R calls (set.seed)
 *
 */

int *Rsample_from(int N, int K, int replace)
{


    /* CHECKS */
    if (K == 0) {
	warning
	    ("Rsample_from.c: K should be > 0; given 0. Behaviour is undefined\n");
	return NULL;
    } else if ((K > N) && !replace) {
	warning
	    ("Rsample_from.c: If K is greater than N, must sample with replacement\n");
	return NULL;
    }

    int *ret = malloc(K * sizeof(double));
    if (!ret) {
	warning
	    ("Rsample_from.c (ERR_MALLOC): not able to allocate double ret of size (%d)",
	     K);
	return NULL;
    }

    int k = 0;
    int totiter = 0;
    GetRNGstate();
    while (k < K) {
	totiter++;
	int p = (int) (unif_rand() * (double) N);
	if (replace) {
	    ret[k] = p;
	    k++;
	} else {
	    int invalid = 0;
	    for (int kk = 0; kk < k; kk++)
		if (p == ret[kk]) {
		    invalid = 1;
		    break;
		}
	    if (!invalid) {
		ret[k] = p;
		k++;
	    }
	}
    }
    PutRNGstate();
    return ret;
}
