#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP frailtyEM_Estep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP frailtyEM_inf_mat_match(SEXP, SEXP, SEXP, SEXP);
extern SEXP frailtyEM_sumxxt(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"frailtyEM_Estep",         (DL_FUNC) &frailtyEM_Estep,         7},
  {"frailtyEM_inf_mat_match", (DL_FUNC) &frailtyEM_inf_mat_match, 4},
  {"frailtyEM_sumxxt",        (DL_FUNC) &frailtyEM_sumxxt,        2},
  {NULL, NULL, 0}
};

void R_init_frailtyEM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
