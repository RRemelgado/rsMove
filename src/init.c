#include R.h
#include Rinternals.h
#include stdlib.h  for NULL
#include R_extRdynload.h

 FIXME 
   Check these declarations against the CFortran source code.


 .Call calls 
extern SEXP _rsMove_intime(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {_rsMove_intime, (DL_FUNC) &_rsMove_intime, 3},
    {NULL, NULL, 0}
};

void R_init_rsMove(DllInfo dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}