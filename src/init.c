#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _IncDTW_BACKTRACK_cpp(SEXP);
extern SEXP _IncDTW_GCM_cpp(SEXP);
extern SEXP _IncDTW_GCM_Sakoe_cpp(SEXP, SEXP);
extern SEXP _IncDTW_IGCM_cpp(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_IGCM_Sakoe_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_IncDTW_BACKTRACK_cpp",  (DL_FUNC) & _IncDTW_BACKTRACK_cpp,  1},
  {"_IncDTW_GCM_cpp",        (DL_FUNC) & _IncDTW_GCM_cpp,        1},
  {"_IncDTW_GCM_Sakoe_cpp",  (DL_FUNC) & _IncDTW_GCM_Sakoe_cpp,  2},
  {"_IncDTW_IGCM_cpp",       (DL_FUNC) & _IncDTW_IGCM_cpp,       3},
  {"_IncDTW_IGCM_Sakoe_cpp", (DL_FUNC) & _IncDTW_IGCM_Sakoe_cpp, 4},
  {NULL, NULL, 0}
};

void R_init_IncDTW(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}