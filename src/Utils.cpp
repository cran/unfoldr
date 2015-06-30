/**
 * utils.cpp
 *
 *  Created on: 01.04.2014
 *      Author: franke
 */

#define BITMAP_GREY "P5"
#define BITMAP_RGB "P6"

#include "Utils.h"
#include <R_ext/Lapack.h>

/* show arguments of .External interface */
SEXP showArgs(SEXP args) {
  int i, nargs;
  Rcomplex cpl;
  const char *name;

  if((nargs = length(args) - 1) > 0) {
    for(i = 0; i < nargs; i++) {
      args = CDR(args);
      name = CHAR(PRINTNAME(TAG(args)));
      switch(TYPEOF(CAR(args))) {
      case REALSXP:
        Rprintf("[%d] '%s' %f\n", i+1, name, REAL(CAR(args))[0]);
        break;
      case LGLSXP:
      case INTSXP:
        Rprintf("[%d] '%s' %d\n", i+1, name, INTEGER(CAR(args))[0]);
        break;
      case CPLXSXP:
        cpl = COMPLEX(CAR(args))[0];
        Rprintf("[%d] '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
        break;
      case STRSXP:
        Rprintf("[%d] '%s' %s\n", i+1, name,
               CHAR(STRING_ELT(CAR(args), 0)));
        break;
      default:
        Rprintf("[%d] '%s' R type\n", i+1, name);
      }
    }
  }
  return(R_NilValue);
}


/* get a single call to an R function */
SEXP getCall(SEXP R_fname, SEXP R_args, SEXP R_rho) {
  SEXP RCallBack = R_NilValue;
  PROTECT(RCallBack = allocVector(LANGSXP, length(R_args)+1 ));
  SETCAR( RCallBack, findFun(install(CHAR(STRING_ELT(R_fname, 0))),R_rho ));

  SEXP p = CDR(RCallBack);
  SEXP names = getAttrib(R_args, R_NamesSymbol);

  for (int i=0; p!=R_NilValue; p=CDR(p),i++) {
    SETCAR(p,VECTOR_ELT(R_args,i));
    SET_TAG(p,install(CHAR(STRING_ELT(names,i))));
  }
  UNPROTECT(1);
  return RCallBack;
}

/* delete the call memebers */
void deleteRCall(R_Calldata call) {
  if(!call)
    return;
  UNPROTECT(call->nprotect);
  Free(call);
}

/* get the elements of a list */
SEXP getListElement (SEXP list, const char *str)
{
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;

}

SEXP getVar(SEXP name, SEXP rho)
{
    SEXP ans;

    if(!isString(name) || length(name) != 1)
        error("name is not a single string");
    if(!isEnvironment(rho))
        error("rho should be an environment");
    ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);
    return ans;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" {

    void real_eval(double *a, int *n, double *evalf, int *err) {
            // size(evalf): n
            // result: evec = a

            int lda = *n,  lwork = 3*lda-1;
            double *work = Calloc(lwork, double);
            F77_NAME(dsyev)("V","U", &lda, a, &lda, evalf, work, &lwork, err);
            Free(work);

    }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Rboolean isNullPtr(SEXP ptr, SEXP type)
{
  if ( TYPEOF(ptr) != EXTPTRSXP ||
      R_ExternalPtrTag(ptr) != type || !R_ExternalPtrAddr(ptr) )
    return TRUE;
  return FALSE;
}

void checkPtr(SEXP ptr, SEXP type)
{
  if ( TYPEOF(ptr) != EXTPTRSXP ||
      R_ExternalPtrTag(ptr) != type || !R_ExternalPtrAddr(ptr) )
    error("Bad Pointer to simulation object");
}

void * getExternalPtr(SEXP ext)
{
  if(!R_ExternalPtrAddr(ext)) {
      warning("Null pointer.\n ");
      return NULL;
  }
  return R_ExternalPtrAddr(ext);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
