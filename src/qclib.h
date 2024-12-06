#ifndef _ECM_H
#define _ECM_H

#include <stdlib.h>

/* Container Triplets representing a clustering solution */
typedef struct Triplets Triplets;
struct Triplets {
    double *prop, *mean, *cov;
};
#define free_Triplets(x) \
	free(x.prop); \
	free(x.mean); \
	free(x.cov);

/* Container for Clustering solution from ECM */
typedef struct ECMout ECMout;
struct ECMout {
    int iter, info;
    int N;
    int K, P;
    int *cluster, *size;
    double eloglik, *tau;
    Triplets params;
};

#define free_ECMout(x) \
	free(x.cluster); \
	free(x.size); \
	free(x.tau); \
	free_Triplets(x.params);


/* Runs the ECM algorithm on data */
ECMout ECM(double *data, int N, int P,
	   int K, double *taus_init, double erc, int itermax,
	   double tol);

/* Runs ECM after initializing weights */
ECMout ECM_winit(double *data, int N, int P,
		 int K, double ecm_erc, int ecm_itermax,
		 double ecm_tol, int kmed_nrestart,
		 int kmed_itermax, double kmed_tol);

/* Kmedian function */
double *kmedians(double *data, int N, int P,
		 int K, int random_nrestart, int itermax,
		 double tol);

/* Random cluster label initializer */
int random_cluster_init(int K);
/* Sample K index out of [0,N] */
int *Rsample_from(int N, int K, int resample);


/* Golden section search (used in ECM) */
void GssERC_v2(double *eigenval, double erc, double eigenval_min,
	       double eigenval_max, double *sumtau, int P,
	       int K);

/* Compute scores using Triplets */
double helper_score_smd(int k, int N, int P,
			double *Delta, double *cov);

double score_smooth(double *data, int N, int P,
		    int K, double *prop, double *mean,
		    double *cov);

double score_hard(double *data, int N, int P,
		  int K, double *prop, double *mean, double *cov);

/* First compute Triplets, then scores */
Triplets cluster_to_triplets(double *data, int N,
			     int P, int K, int *cluster);
double score_smooth_cluster(double *data, int N,
			    int P, int K, int *cluster);
double score_hard_cluster(double *data, int N,
			  int P, int K, int *cluster);

//double        score_hard_cluster(double *data, int N, int P, int K, double *cluster);


/* Quick select algorithm */
double quickselect(double *vector_extract_k, int vec_size,
		   int extract_this_element);


#endif				/* Include only if not included already _ECM_H */
