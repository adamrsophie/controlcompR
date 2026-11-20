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
int cc_square_real(SEXP x, const char *name);
int cc_mode(SEXP mode_);
double cc_tol(SEXP tol_);
int cc_lwork(double wq, int info, const char *routine);

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

/* want_ = c(Q, Z): which Schur vector matrices to compute; a matrix that
 * is not wanted is returned as NULL and LAPACK does not form it. */
SEXP cc_qz(SEXP a_, SEXP b_, SEXP mode_, SEXP tol_, SEXP want_)
{
    if (!isLogical(want_) || length(want_) != 2 ||
        LOGICAL(want_)[0] == NA_LOGICAL || LOGICAL(want_)[1] == NA_LOGICAL)
        error("`want` must be two non-missing logical values");
    int wq_ = LOGICAL(want_)[0], wz_ = LOGICAL(want_)[1];
    int n = cc_square_real(a_, "A"), mode = cc_mode(mode_), info = 0,
        sdim = 0, lwork = -1, ld = n > 0 ? n : 1;
    if (cc_square_real(b_, "B") != n) error("`A` and `B` must have the same size");
    cc_sel_tol = cc_tol(tol_);
    SEXP AA = PROTECT(duplicate(a_));
    SEXP BB = PROTECT(duplicate(b_));
    SEXP Q = PROTECT(wq_ ? allocMatrix(REALSXP, n, n) : allocVector(REALSXP, 1));
    SEXP Z = PROTECT(wz_ ? allocMatrix(REALSXP, n, n) : allocVector(REALSXP, 1));
    const char *jl = wq_ ? "V" : "N", *jr = wz_ ? "V" : "N";
    int ldq = wq_ ? ld : 1, ldz = wz_ ? ld : 1;
    SEXP ar = PROTECT(allocVector(REALSXP, n));
    SEXP ai = PROTECT(allocVector(REALSXP, n));
    SEXP be = PROTECT(allocVector(REALSXP, n));
    int *bwork = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
    const char *sort = mode == 0 ? "N" : "S";
    int (*sel)(double *, double *, double *) = mode == 2 ? gsel_dt : gsel_ct;
    double wq = 0.0;
    F77_CALL(dgges)(jl, jr, sort, sel, &n, REAL(AA), &ld, REAL(BB), &ld, &sdim,
                    REAL(ar), REAL(ai), REAL(be), REAL(Q), &ldq, REAL(Z), &ldz,
                    &wq, &lwork, bwork, &info CC_FCONE CC_FCONE CC_FCONE);
    lwork = cc_lwork(wq, info, "dgges");
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgges)(jl, jr, sort, sel, &n, REAL(AA), &ld, REAL(BB), &ld, &sdim,
                    REAL(ar), REAL(ai), REAL(be), REAL(Q), &ldq, REAL(Z), &ldz,
                    work, &lwork, bwork, &info CC_FCONE CC_FCONE CC_FCONE);
    if (info < 0) error("LAPACK dgges: argument %d is invalid", -info);
    const char *nms[] = {"AA", "BB", "Q", "Z", "alphar", "alphai", "beta", "sdim", "info"};
    SEXP out = PROTECT(cc_named_list(9, nms));
    SET_VECTOR_ELT(out, 0, AA);
    SET_VECTOR_ELT(out, 1, BB);
    SET_VECTOR_ELT(out, 2, wq_ ? Q : R_NilValue);
    SET_VECTOR_ELT(out, 3, wz_ ? Z : R_NilValue);
    SET_VECTOR_ELT(out, 4, ar);
    SET_VECTOR_ELT(out, 5, ai);
    SET_VECTOR_ELT(out, 6, be);
    SET_VECTOR_ELT(out, 7, ScalarInteger(sdim));
    SET_VECTOR_ELT(out, 8, ScalarInteger(info));
    UNPROTECT(8);
    return out;
}

