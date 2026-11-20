/* Solve T X T' - X + Q = 0 for quasi upper triangular T (real Schur form),
 * column block by column block from the right, with the same block
 * recursion as the reference R implementation. Blocks are given by their
 * 0-based starts and sizes (1 or 2). */
#include <R.h>
#include <Rinternals.h>
#include <math.h>

int cc_square_real(SEXP x, const char *name);

/* Solve the k x k system K x = b (k <= 4) by Gaussian elimination with
 * partial pivoting; K and b are overwritten. Returns 0 or 1 if singular. */
static int small_solve(int k, double *K, double *b)
{
    for (int c = 0; c < k; c++) {
        int piv = c;
        for (int r = c + 1; r < k; r++)
            if (fabs(K[r + 4 * c]) > fabs(K[piv + 4 * c])) piv = r;
        if (K[piv + 4 * c] == 0.0) return 1;
        if (piv != c) {
            for (int j = 0; j < k; j++) {
                double t = K[c + 4 * j]; K[c + 4 * j] = K[piv + 4 * j]; K[piv + 4 * j] = t;
            }
            double t = b[c]; b[c] = b[piv]; b[piv] = t;
        }
        for (int r = c + 1; r < k; r++) {
            double f = K[r + 4 * c] / K[c + 4 * c];
            for (int j = c; j < k; j++) K[r + 4 * j] -= f * K[c + 4 * j];
            b[r] -= f * b[c];
        }
    }
    for (int c = k - 1; c >= 0; c--) {
        double s = b[c];
        for (int j = c + 1; j < k; j++) s -= K[c + 4 * j] * b[j];
        b[c] = s / K[c + 4 * c];
    }
    return 0;
}

SEXP cc_stein(SEXP t_, SEXP q_, SEXP starts_, SEXP sizes_)
{
    int n = cc_square_real(t_, "T");
    if (cc_square_real(q_, "Q") != n) error("`T` and `Q` must have the same size");
    if (!isInteger(starts_) || !isInteger(sizes_) || length(starts_) != length(sizes_))
        error("`starts` and `sizes` must be integer vectors of equal length");
    int p = length(starts_);
    const int *bs = INTEGER(starts_), *bz = INTEGER(sizes_);
    int covered = 0;
    for (int b = 0; b < p; b++) {
        if (bs[b] != covered || (bz[b] != 1 && bz[b] != 2))
            error("invalid Schur block structure");
        covered += bz[b];
    }
    if (covered != n) error("invalid Schur block structure");
    const double *T = REAL(t_), *Q = REAL(q_);
    SEXP x_ = PROTECT(allocMatrix(REALSXP, n, n));
    double *X = REAL(x_);
    for (R_xlen_t i = 0; i < (R_xlen_t) n * n; i++) X[i] = 0.0;
    double *Y = (double *) R_alloc(n > 0 ? 2 * (size_t) n : 1, sizeof(double));
    double *W = (double *) R_alloc(n > 0 ? 2 * (size_t) n : 1, sizeof(double));
#define TT(i, j) T[(i) + (R_xlen_t) n * (j)]
#define XX(i, j) X[(i) + (R_xlen_t) n * (j)]
    for (int jb = p - 1; jb >= 0; jb--) {
        int j0 = bs[jb], nj = bz[jb], later = j0 + nj;
        /* Y = X[, later:] T[J, later:]' */
        for (int c = 0; c < nj; c++)
            for (int r = 0; r < n; r++) {
                double s = 0.0;
                for (int l = later; l < n; l++) s += XX(r, l) * TT(j0 + c, l);
                Y[r + (R_xlen_t) n * c] = s;
            }
        for (int ib = p - 1; ib >= 0; ib--) {
            int i0 = bs[ib], ni = bz[ib], gt = i0 + ni, k = ni * nj;
            /* W[r, c] = (X[r, J] Tjj')[c] for the already solved rows r >= gt */
            for (int c = 0; c < nj; c++)
                for (int r = gt; r < n; r++) {
                    double s = 0.0;
                    for (int d = 0; d < nj; d++) s += XX(r, j0 + d) * TT(j0 + c, j0 + d);
                    W[r + (R_xlen_t) n * c] = s;
                }
            double rhs[4], K[16];
            for (int c = 0; c < nj; c++)
                for (int a = 0; a < ni; a++) {
                    double s = -Q[(i0 + a) + (R_xlen_t) n * (j0 + c)];
                    for (int r = i0; r < n; r++) s -= TT(i0 + a, r) * Y[r + (R_xlen_t) n * c];
                    for (int r = gt; r < n; r++) s -= TT(i0 + a, r) * W[r + (R_xlen_t) n * c];
                    rhs[a + ni * c] = s;
                }
            /* (Tjj kron Tii - I) vec(Xij) = vec(rhs) */
            for (int c = 0; c < nj; c++)
                for (int a = 0; a < ni; a++)
                    for (int d = 0; d < nj; d++)
                        for (int b = 0; b < ni; b++)
                            K[(a + ni * c) + 4 * (b + ni * d)] =
                                TT(j0 + c, j0 + d) * TT(i0 + a, i0 + b) -
                                ((a == b && c == d) ? 1.0 : 0.0);
            if (small_solve(k, K, rhs)) {
                UNPROTECT(1);
                error("the Stein equation is singular (eigenvalue products equal to 1)");
            }
            for (int c = 0; c < nj; c++)
                for (int a = 0; a < ni; a++) XX(i0 + a, j0 + c) = rhs[a + ni * c];
        }
    }
#undef TT
#undef XX
    UNPROTECT(1);
    return x_;
}
