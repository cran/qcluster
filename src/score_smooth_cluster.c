#include "qclib.h"
/*
 *  score_smooth_cluster
 *
 *  Computes QS with given data and cluster labels
 *
 *
 *  Arguments:
 *  - [in] double *data : pointer to data (a (P,N) matrix)
 *  - [in] int N        : number of observations in data
 *  - [in] int P        : number of dimensions (features) in data
 *  - [in] int K        : number of clusters in solution
 *  - [in] int *cluster : cluster labels
 *
 *  Returns double:
 *  the smooth scoring computed on data, using the clustering solution
 *  given by "cluster"
 *
 */
double score_smooth_cluster(double *data, int N, int P, int K,
			    int *cluster)
{

    Triplets params = cluster_to_triplets(data, N, P, K, cluster);
    double res =
	score_smooth(data, N, P, K, params.prop, params.mean, params.cov);

free_Triplets(params);
    return res;
}
