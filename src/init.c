#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _IncDTW_BACKTRACK_cpp(SEXP);
extern SEXP _IncDTW_BACKTRACK2II_cpp(SEXP, SEXP);
extern SEXP _IncDTW_BACKTRACK2IN_cpp(SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec(SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ea(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_inc_ws(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_v32(SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ws(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_cpp_dtw2vec_ws_ea(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IncDTW_GCM_cpp(SEXP);
extern SEXP _IncDTW_GCM_Sakoe_cpp(SEXP, SEXP);
extern SEXP _IncDTW_IGCM_cpp(SEXP, SEXP, SEXP);
extern SEXP _IncDTW_IGCM_Sakoe_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
   {"_IncDTW_BACKTRACK_cpp",      (DL_FUNC) &_IncDTW_BACKTRACK_cpp,      1},
   {"_IncDTW_BACKTRACK2II_cpp",   (DL_FUNC) &_IncDTW_BACKTRACK2II_cpp,   2},
   {"_IncDTW_BACKTRACK2IN_cpp",   (DL_FUNC) &_IncDTW_BACKTRACK2IN_cpp,   2},
   {"_IncDTW_cpp_dtw2vec",        (DL_FUNC) &_IncDTW_cpp_dtw2vec,        2},
   {"_IncDTW_cpp_dtw2vec_ea",     (DL_FUNC) &_IncDTW_cpp_dtw2vec_ea,     3},
   {"_IncDTW_cpp_dtw2vec_inc",    (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc,    3},
   {"_IncDTW_cpp_dtw2vec_inc_ws", (DL_FUNC) &_IncDTW_cpp_dtw2vec_inc_ws, 5},
   {"_IncDTW_cpp_dtw2vec_v32",    (DL_FUNC) &_IncDTW_cpp_dtw2vec_v32,    2},
   {"_IncDTW_cpp_dtw2vec_ws",     (DL_FUNC) &_IncDTW_cpp_dtw2vec_ws,     3},
   {"_IncDTW_cpp_dtw2vec_ws_ea",  (DL_FUNC) &_IncDTW_cpp_dtw2vec_ws_ea,  4},
   {"_IncDTW_GCM_cpp",            (DL_FUNC) &_IncDTW_GCM_cpp,            1},
   {"_IncDTW_GCM_Sakoe_cpp",      (DL_FUNC) &_IncDTW_GCM_Sakoe_cpp,      2},
   {"_IncDTW_IGCM_cpp",           (DL_FUNC) &_IncDTW_IGCM_cpp,           3},
   {"_IncDTW_IGCM_Sakoe_cpp",     (DL_FUNC) &_IncDTW_IGCM_Sakoe_cpp,     4},
   {NULL, NULL, 0}
};

void R_init_IncDTW(DllInfo *dll)
{
   R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
}