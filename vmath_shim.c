/* Portable vector-math shim so amica15 builds without a vendor math library.
 *
 * amica15.f90's non-MKL branch calls AMD LibM's vectorized exp/log
 * (vrda_exp / vrda_log); its MKL branch uses Intel VML (vdExp / vdLn). Neither
 * vendor library is present in a plain gfortran + LAPACK build (e.g. on Apple
 * Silicon, or a generic Linux box without ACML), so the link fails on
 * vrda_exp_ / vrda_log_.
 *
 * These are simple element-wise ops, so this shim provides them as loops over
 * standard libm exp/log. gfortran calls them by reference with a trailing
 * underscore (`call vrda_exp(n, x, y)` -> `vrda_exp_(int* n, double* x,
 * double* y)`), which is the ABI matched here. Scalar libm exp/log are
 * IEEE-accurate; a vendor SIMD math library would only be faster, not more
 * accurate.
 */
#include <math.h>

void vrda_exp_(const int *n, const double *x, double *y) {
    int i, m = *n;
    for (i = 0; i < m; i++) y[i] = exp(x[i]);
}

void vrda_log_(const int *n, const double *x, double *y) {
    int i, m = *n;
    for (i = 0; i < m; i++) y[i] = log(x[i]);
}
