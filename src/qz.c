/* Ordered generalized real Schur decomposition (QZ) through LAPACK dgges.
 * R's R_ext/Lapack.h declares dgges without its SDIM argument, so this file
 * declares the reference LAPACK prototype itself and does not include that
 * header. */
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <math.h>
#ifdef FC_LEN_T
# define CC_FCLEN , FC_LEN_T
# define CC_FCONE , (FC_LEN_T) 1
#else
# define CC_FCLEN
# define CC_FCONE
#endif

extern double cc_sel_tol;
SEXP cc_named_list(int n, const char **names);

extern void F77_NAME(dgges)(const char *jobvsl, const char *jobvsr,
                            const char *sort,
                            int (*selctg)(double *, double *, double *),
                            const int *n, double *a, const int *lda,
                            double *b, const int *ldb, int *sdim,
                            double *alphar, double *alphai, double *beta,
                            double *vsl, const int *ldvsl, double *vsr,
                            const int *ldvsr, double *work, const int *lwork,
                            int *bwork, int *info CC_FCLEN CC_FCLEN CC_FCLEN);

static int gsel_ct(double *ar, double *ai, double *b)
{
    (void) ai;
    /* lambda = (ar + i ai) / b; infinite eigenvalues (b == 0) are not stable */
    if (*b == 0.0) return 0;
    return (*ar / *b) < -cc_sel_tol;
}

static int gsel_dt(double *ar, double *ai, double *b)
{
    return hypot(*ar, *ai) < (1.0 - cc_sel_tol) * fabs(*b);
}

SEXP cc_qz(SEXP a_, SEXP b_, SEXP mode_, SEXP tol_)
{
    int n = nrows(a_), mode = asInteger(mode_), info = 0, sdim = 0, lwork = -1;
    cc_sel_tol = asReal(tol_);
    SEXP AA = PROTECT(duplicate(a_));
    SEXP BB = PROTECT(duplicate(b_));
    SEXP Q = PROTECT(allocMatrix(REALSXP, n, n));
    SEXP Z = PROTECT(allocMatrix(REALSXP, n, n));
    SEXP ar = PROTECT(allocVector(REALSXP, n));
    SEXP ai = PROTECT(allocVector(REALSXP, n));
    SEXP be = PROTECT(allocVector(REALSXP, n));
    int *bwork = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
    const char *sort = mode == 0 ? "N" : "S";
    int (*sel)(double *, double *, double *) = mode == 2 ? gsel_dt : gsel_ct;
    double wq;
    F77_CALL(dgges)("V", "V", sort, sel, &n, REAL(AA), &n, REAL(BB), &n, &sdim,
                    REAL(ar), REAL(ai), REAL(be), REAL(Q), &n, REAL(Z), &n,
                    &wq, &lwork, bwork, &info CC_FCONE CC_FCONE CC_FCONE);
    lwork = (int) wq;
    if (lwork < 1) lwork = 1;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgges)("V", "V", sort, sel, &n, REAL(AA), &n, REAL(BB), &n, &sdim,
                    REAL(ar), REAL(ai), REAL(be), REAL(Q), &n, REAL(Z), &n,
                    work, &lwork, bwork, &info CC_FCONE CC_FCONE CC_FCONE);
    const char *nms[] = {"AA", "BB", "Q", "Z", "alphar", "alphai", "beta", "sdim", "info"};
    SEXP out = PROTECT(cc_named_list(9, nms));
    SET_VECTOR_ELT(out, 0, AA);
    SET_VECTOR_ELT(out, 1, BB);
    SET_VECTOR_ELT(out, 2, Q);
    SET_VECTOR_ELT(out, 3, Z);
    SET_VECTOR_ELT(out, 4, ar);
    SET_VECTOR_ELT(out, 5, ai);
    SET_VECTOR_ELT(out, 6, be);
    SET_VECTOR_ELT(out, 7, ScalarInteger(sdim));
    SET_VECTOR_ELT(out, 8, ScalarInteger(info));
    UNPROTECT(8);
    return out;
}

