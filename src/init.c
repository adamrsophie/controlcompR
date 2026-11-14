#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP cc_schur(SEXP, SEXP, SEXP);
SEXP cc_qz(SEXP, SEXP, SEXP, SEXP);
SEXP cc_trsyl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cc_schur", (DL_FUNC) &cc_schur, 3},
    {"cc_qz", (DL_FUNC) &cc_qz, 4},
    {"cc_trsyl", (DL_FUNC) &cc_trsyl, 6},
    {NULL, NULL, 0}
};

void R_init_controlcompR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
