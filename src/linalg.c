/* Thin wrappers around LAPACK routines that base R does not expose:
 * ordered real Schur (dgees), ordered generalized Schur (dgges) and the
 * quasi-triangular Sylvester solver (dtrsyl). */
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <math.h>
#ifndef FCONE
# define FCONE
#endif

/* Selection criteria for the leading block of the ordered Schur form.
 * Mode 1: continuous-time stable (real part < -tol).
 * Mode 2: discrete-time stable (modulus < 1 - tol). */
double cc_sel_tol = 0.0;

static int sel_ct(const double *wr, const double *wi)
{
    (void) wi;
    return *wr < -cc_sel_tol;
}

static int sel_dt(const double *wr, const double *wi)
{
    return hypot(*wr, *wi) < 1.0 - cc_sel_tol;
}

SEXP cc_named_list(int n, const char **names)
{
    SEXP out = PROTECT(allocVector(VECSXP, n));
    SEXP nms = PROTECT(allocVector(STRSXP, n));
    for (int i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
}

SEXP cc_schur(SEXP a_, SEXP mode_, SEXP tol_)
{
    int n = nrows(a_), mode = asInteger(mode_), info = 0, sdim = 0, lwork = -1;
    cc_sel_tol = asReal(tol_);
    SEXP T = PROTECT(duplicate(a_));
    SEXP Z = PROTECT(allocMatrix(REALSXP, n, n));
    SEXP wr = PROTECT(allocVector(REALSXP, n));
    SEXP wi = PROTECT(allocVector(REALSXP, n));
    int *bwork = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
    const char *sort = mode == 0 ? "N" : "S";
    int (*sel)(const double *, const double *) = mode == 2 ? sel_dt : sel_ct;
    double wq;
    F77_CALL(dgees)("V", sort, sel, &n, REAL(T), &n, &sdim, REAL(wr), REAL(wi),
                    REAL(Z), &n, &wq, &lwork, bwork, &info FCONE FCONE);
    lwork = (int) wq;
    if (lwork < 1) lwork = 1;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgees)("V", sort, sel, &n, REAL(T), &n, &sdim, REAL(wr), REAL(wi),
                    REAL(Z), &n, work, &lwork, bwork, &info FCONE FCONE);
    const char *nms[] = {"T", "Z", "wr", "wi", "sdim", "info"};
    SEXP out = PROTECT(cc_named_list(6, nms));
    SET_VECTOR_ELT(out, 0, T);
    SET_VECTOR_ELT(out, 1, Z);
    SET_VECTOR_ELT(out, 2, wr);
    SET_VECTOR_ELT(out, 3, wi);
    SET_VECTOR_ELT(out, 4, ScalarInteger(sdim));
    SET_VECTOR_ELT(out, 5, ScalarInteger(info));
    UNPROTECT(5);
    return out;
}

/* Solve op(A) X + isgn * X op(B) = scale * C with A, B quasi-triangular. */
SEXP cc_trsyl(SEXP a_, SEXP b_, SEXP c_, SEXP trana_, SEXP tranb_, SEXP isgn_)
{
    int m = nrows(a_), n = nrows(b_), isgn = asInteger(isgn_), info = 0;
    double scale = 1.0;
    const char *ta = CHAR(STRING_ELT(trana_, 0));
    const char *tb = CHAR(STRING_ELT(tranb_, 0));
    SEXP X = PROTECT(duplicate(c_));
    F77_CALL(dtrsyl)(ta, tb, &isgn, &m, &n, REAL(a_), &m, REAL(b_), &n,
                     REAL(X), &m, &scale, &info FCONE FCONE);
    const char *nms[] = {"X", "scale", "info"};
    SEXP out = PROTECT(cc_named_list(3, nms));
    SET_VECTOR_ELT(out, 0, X);
    SET_VECTOR_ELT(out, 1, ScalarReal(scale));
    SET_VECTOR_ELT(out, 2, ScalarInteger(info));
    UNPROTECT(2);
    return out;
}
