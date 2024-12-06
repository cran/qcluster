/* Golden Section Search for ERC */

#include <float.h>
#include <math.h>
#include <stdio.h>

#define GRCINV (pow(5.0, 0.5) - 1)/2
#define TOL pow(DBL_EPSILON, 0.5)
#define ITERMAX 99

/**
 *
 * name: func
 * @param[in] x domain value
 * @param[in] eigenval pointer to matrix of eigenvalues
 * @param[in] sumtau pointer to array of posterior sums
 * @param[in] P the number of dimensions
 * @param[in] K the number of cluster
 * @return f the value of the function evaluated at x
 *
 */
static double
func(double x, double erc, double *eigenval, double *sumtau, int P, int K)
{

    double res = 0;

    x = exp(x);
    for (int k = 0; k < K; k++) {
        double tmp = 0;
        for (int p = 0; p < P; p++) {
            if (*eigenval < x) {
                tmp += log(x) + *eigenval / x;
            } else if (*eigenval > x * erc) {
                tmp += log(x * erc) + *eigenval / (x * erc);
            } else {
                tmp += log(*eigenval) + 1;
            }
            eigenval++;
        }
        res += tmp * *sumtau++;
    }

    return res;
}

/**
 * Golden section search for the ERC
 * name: GssERC_v2
 * @param[in,out] eigenval the (K,L) matrix of eigenvalues
 * @param[in] erc eigenratio constraint
 * @param[] eigenval_min minimum value of eigenval
 * @param[] eigenval_max maximum value of eigenval
 * @param[] sumtau the K vector of row sums of tau (posterior) weights
 * @param[in] P the number of dimensions
 * @param[in] K the number of cluster
 */
void
GssERC_v2(double *eigenval, double erc, double eigenval_min,
          double eigenval_max, double *sumtau, int P, int K)
{

    int iter = 0;
    double a, b, c, d, fc, fd;
    a = eigenval_min > DBL_EPSILON ? log(eigenval_min) : log(DBL_EPSILON);
    b = eigenval_max < DBL_MAX ? log(eigenval_max) : log(DBL_MAX);

    /* Initial bracketing */
    d = a + GRCINV * (b - a);
    c = b - GRCINV * (b - a);

    /* Eval f(c) */
    iter++;
    fc = func(c, erc, eigenval, sumtau, P, K);

    /* Eval f(d) */
    iter++;
    fd = func(d, erc, eigenval, sumtau, P, K);

    while ((fabs(b - a) > TOL) && (iter++ < ITERMAX)) {
        if (fc > fd) {
            a = c;
            c = d;
            d = a + GRCINV * (b - a);
            fc = fd;
            fd = func(d, erc, eigenval, sumtau, P, K);
        } else {
            b = d;
            d = c;
            c = b - GRCINV * (b - a);
            fd = fc;
            fc = func(c, erc, eigenval, sumtau, P, K);
        }
    }

    // we now use d as a temporary container
    d = exp((a + b) / 2);
    for (int k = 0; k < K; k++) {
        for (int p = 0; p < P; p++) {
            if (*eigenval < d)
                *eigenval = d;
            else if (*eigenval > d * erc)
                *eigenval = d * erc;
            eigenval++;
        }
    }
}
