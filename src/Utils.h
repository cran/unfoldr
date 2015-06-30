#ifndef UTILS_H_
#define UTILS_H_

#include "Rheaders.h"

/* functions with R types */

void *getExternalPtr(SEXP ext);
void checkPtr(SEXP ptr,SEXP type);
Rboolean isNullPtr(SEXP ptr, SEXP type);

SEXP getVar(SEXP name, SEXP rho);
SEXP getListElement (SEXP list, const char *str);

/* internally used */

#ifdef __cplusplus
extern "C" {
#endif

/*  real eigenvalue decomposition */
void real_eval(double *a, int *n, double *evalf, int *err);
/* pointer to some R's random generators */
typedef double (*rdist2_t)(double, double);
/* a const dummy function */
inline double rconst(double x, double dummy=0) { return x; }


/*
 * R call data struct, fname,args and call could also be
 * lists of names, args and calls to functions
 */
typedef struct R_Calldata_s {
    SEXP fname,args,rho,call;
    int nprotect;
} *R_Calldata;


#ifdef __cplusplus
}
#endif

template<typename R_TYPE>
struct R_eval_t {
  typedef R_TYPE return_type;
  SEXP call, rho;
  R_eval_t(SEXP _call, SEXP _rho) :  call(_call), rho(_rho)  {};
  inline return_type operator()() { return AS_NUMERIC(eval(call,rho));  }
};

template<>
struct R_eval_t<double> {
  SEXP call, rho;
  R_eval_t(SEXP _call, SEXP _rho) :  call(_call), rho(_rho)
  {};
  inline double operator()() {
    return asReal(eval(call,rho));
  }
};

template<typename F >
struct R_rndGen_t {
  double p,q;
  F fn;
  R_rndGen_t(double _p,double _q, F _fn) : p(_p), q(_q), fn(_fn) {};
  inline double operator()() { return fn(p,q); }
};

/**
 * \brief Show aruments of R function call
 *
 * @param args list of arguments
 * @return
 */
SEXP showArgs(SEXP args);

/**
 * \brief Define R call to user function from C level
 *
 * @param  R_fname (string) name of function to call
 * @param  R_args  list of arguments function
 * @param  R_rho   R environment
 * @return SEXP call as evaluated by eval(Call, Env)
 */
SEXP getCall(SEXP R_fname, SEXP R_args, SEXP R_rho);

/**
 * \brief Delete R_Calldata struct and and UNPROTECT SEXPs
 * @param call
 */
void deleteRCall(R_Calldata call);


#endif /* UTILS_H_ */
