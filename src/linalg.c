/* Thin wrappers around LAPACK routines that base R does not expose:
 * ordered real Schur (dgees), ordered generalized Schur (dgges) and the
 * quasi-triangular Sylvester solver (dtrsyl). */
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <math.h>
#include <limits.h>
#include <string.h>
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

/* Entry-point checks. The R wrappers already pass valid arguments; these
 * checks make a direct .Call with the wrong types or dimensions fail with
 * an R error instead of reading or writing out of bounds. */
int cc_square_real(SEXP x, const char *name)
{
    if (!isReal(x) || !isMatrix(x))
        error("`%s` must be a double matrix", name);
    int n = nrows(x);
    if (ncols(x) != n) error("`%s` must be square", name);
    return n;
}

int cc_mode(SEXP mode_)
{
    if (length(mode_) != 1) error("`mode` must be a single integer");
    int mode = asInteger(mode_);
    if (mode == NA_INTEGER || mode < 0 || mode > 2)
        error("`mode` must be 0, 1 or 2");
    return mode;
}

double cc_tol(SEXP tol_)
{
    if (length(tol_) != 1) error("`tol` must be a single number");
    double tol = asReal(tol_);
    if (!R_FINITE(tol)) error("`tol` must be finite");
    return tol;
}

/* Workspace size from a LAPACK query, checked before the cast to int. */
int cc_lwork(double wq, int info, const char *routine)
{
    if (info != 0) error("LAPACK %s workspace query failed (info %d)", routine, info);
    if (!R_FINITE(wq) || wq > (double) INT_MAX)
        error("LAPACK %s workspace size is not representable", routine);
    int lwork = (int) wq;
    return lwork < 1 ? 1 : lwork;
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
    int n = cc_square_real(a_, "A"), mode = cc_mode(mode_), info = 0,
        sdim = 0, lwork = -1, ld = n > 0 ? n : 1;
    cc_sel_tol = cc_tol(tol_);
    SEXP T = PROTECT(duplicate(a_));
    SEXP Z = PROTECT(allocMatrix(REALSXP, n, n));
    SEXP wr = PROTECT(allocVector(REALSXP, n));
    SEXP wi = PROTECT(allocVector(REALSXP, n));
    int *bwork = (int *) R_alloc(n > 0 ? n : 1, sizeof(int));
    const char *sort = mode == 0 ? "N" : "S";
    int (*sel)(const double *, const double *) = mode == 2 ? sel_dt : sel_ct;
    double wq = 0.0;
    F77_CALL(dgees)("V", sort, sel, &n, REAL(T), &ld, &sdim, REAL(wr), REAL(wi),
                    REAL(Z), &ld, &wq, &lwork, bwork, &info FCONE FCONE);
    lwork = cc_lwork(wq, info, "dgees");
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgees)("V", sort, sel, &n, REAL(T), &ld, &sdim, REAL(wr), REAL(wi),
                    REAL(Z), &ld, work, &lwork, bwork, &info FCONE FCONE);
    if (info < 0) error("LAPACK dgees: argument %d is invalid", -info);
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
    int m = cc_square_real(a_, "A"), n = cc_square_real(b_, "B"), info = 0;
    if (!isReal(c_) || !isMatrix(c_) || nrows(c_) != m || ncols(c_) != n)
        error("`C` must be a double matrix with nrow(A) rows and nrow(B) columns");
    if (!isString(trana_) || length(trana_) != 1 ||
        !isString(tranb_) || length(tranb_) != 1)
        error("`trana` and `tranb` must be single strings");
    const char *ta = CHAR(STRING_ELT(trana_, 0));
    const char *tb = CHAR(STRING_ELT(tranb_, 0));
    if ((strcmp(ta, "N") && strcmp(ta, "T")) || (strcmp(tb, "N") && strcmp(tb, "T")))
        error("`trana` and `tranb` must be \"N\" or \"T\"");
    if (length(isgn_) != 1) error("`isgn` must be 1 or -1");
    int isgn = asInteger(isgn_);
    if (isgn != 1 && isgn != -1) error("`isgn` must be 1 or -1");
    int lda = m > 0 ? m : 1, ldb = n > 0 ? n : 1;
    double scale = 1.0;
    SEXP X = PROTECT(duplicate(c_));
    F77_CALL(dtrsyl)(ta, tb, &isgn, &m, &n, REAL(a_), &lda, REAL(b_), &ldb,
                     REAL(X), &lda, &scale, &info FCONE FCONE);
    if (info < 0) error("LAPACK dtrsyl: argument %d is invalid", -info);
    const char *nms[] = {"X", "scale", "info"};
    SEXP out = PROTECT(cc_named_list(3, nms));
    SET_VECTOR_ELT(out, 0, X);
    SET_VECTOR_ELT(out, 1, ScalarReal(scale));
    SET_VECTOR_ELT(out, 2, ScalarInteger(info));
    UNPROTECT(2);
    return out;
}
