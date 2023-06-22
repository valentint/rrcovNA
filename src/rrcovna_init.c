#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(emnint)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fnanmcd)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
                              void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
                              void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
                              void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
                              void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"emnint",  (DL_FUNC) &F77_NAME(emnint),  11},
    {"fnanmcd", (DL_FUNC) &F77_NAME(fnanmcd), 42},
    {NULL, NULL, 0}
};

void R_init_rrcovNA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

