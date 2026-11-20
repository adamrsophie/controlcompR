/* Frequency response from a complex Schur form: for each point z,
 * Y = Ct (z I - T)^-1 Bt with T upper triangular, by back substitution
 * (O(n^2 m) per point after one O(n^3) factorization). A point within the
 * relative distance `tol` of an eigenvalue (a pole to rounding accuracy)
 * or where the result is not finite gives NA. */
#include <R.h>
#include <Rinternals.h>
#include <complex.h>
#include <math.h>
#include <limits.h>

static void check_cplx(SEXP x, int nr, int nc, const char *name)
{
    if (!isComplex(x) || !isMatrix(x) || nrows(x) != nr || ncols(x) != nc)
        error("`%s` must be a %d x %d complex matrix", name, nr, nc);
}

SEXP cc_resolvent(SEXP t_, SEXP b_, SEXP c_, SEXP z_, SEXP tol_)
{
    if (!isReal(tol_) || length(tol_) != 1 || !R_FINITE(REAL(tol_)[0]) || REAL(tol_)[0] < 0)
        error("`tol` must be a nonnegative finite number");
    double tol = REAL(tol_)[0];
    if (!isComplex(t_) || !isMatrix(t_) || nrows(t_) != ncols(t_))
        error("`T` must be a square complex matrix");
    int n = nrows(t_);
    if (!isComplex(b_) || !isMatrix(b_)) error("`B` must be a complex matrix");
    int m = ncols(b_);
    check_cplx(b_, n, m, "B");
    if (!isComplex(c_) || !isMatrix(c_)) error("`C` must be a complex matrix");
    int p = nrows(c_);
    check_cplx(c_, p, n, "C");
    if (!isComplex(z_)) error("`z` must be a complex vector");
    R_xlen_t nz = XLENGTH(z_);
    const double complex *T = (const double complex *) COMPLEX(t_);
    const double complex *B = (const double complex *) COMPLEX(b_);
    const double complex *C = (const double complex *) COMPLEX(c_);
    const double complex *Z = (const double complex *) COMPLEX(z_);
    SEXP dims = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dims)[0] = p; INTEGER(dims)[1] = m; INTEGER(dims)[2] = (int) nz;
    if ((double) nz > INT_MAX) error("too many frequency points");
    SEXP out = PROTECT(allocArray(CPLXSXP, dims));
    double complex *Y = (double complex *) COMPLEX(out);
    double complex *X = (double complex *) R_alloc(n > 0 ? (size_t) n : 1,
                                                   sizeof(double complex));
    for (R_xlen_t k = 0; k < nz; k++) {
        double complex z = Z[k];
        int bad = 0;
        for (int j = 0; j < m && !bad; j++) {
            for (int i = n - 1; i >= 0; i--) {
                double complex s = B[i + (R_xlen_t) n * j];
                for (int l = i + 1; l < n; l++) s += T[i + (R_xlen_t) n * l] * X[l];
                double complex d = z - T[i + (R_xlen_t) n * i];
                if (cabs(d) <= tol * fmax(cabs(z), cabs(T[i + (R_xlen_t) n * i]))) {
                    bad = 1; break;
                }
                X[i] = s / d;
            }
            if (bad) break;
            for (int r = 0; r < p; r++) {
                double complex s = 0;
                for (int l = 0; l < n; l++) s += C[r + (R_xlen_t) p * l] * X[l];
                if (!isfinite(creal(s)) || !isfinite(cimag(s))) { bad = 1; break; }
                Y[r + (R_xlen_t) p * j + (R_xlen_t) p * m * k] = s;
            }
        }
        if (bad) {
            for (int j = 0; j < m; j++)
                for (int r = 0; r < p; r++) {
                    Rcomplex na; na.r = NA_REAL; na.i = NA_REAL;
                    COMPLEX(out)[r + (R_xlen_t) p * j + (R_xlen_t) p * m * k] = na;
                }
        }
    }
    UNPROTECT(2);
    return out;
}
