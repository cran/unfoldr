/**
 * utils.h
 *
 */

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

/* pointer to some R's random generators */
typedef double (*rdist2_t)(double, double);
/* a const dummy function */
inline double rconst(double x, double dummy=0) { return x; }


/*
 * R call data struct, fname,args and call could also be
 * lists of names, args and calls to functions
 */
typedef struct R_Calldata_s {
    SEXP fname,args,rho,label,call;
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
  {
  };
  inline double operator()() {
    return asReal(eval(call,rho));
  }
};

//template<typename F >
//struct R_rndGen_t {
//  double mx,sdx;
//  F fn;
//  R_rndGen_t(double p,double q, F rdist) : mx(p), sdx(q), fn(rdist) {};
//  inline double operator()() { return fn(mx,sdx); }
//};


struct R_rlnorm_t {
  double mx,sdx;
  R_rlnorm_t(double p,double q)
    : mx(p), sdx(q)
  {};

  inline double operator()() {  return rlnorm(mx,sdx); }
};

template<typename F >
struct R_rndGen_t {
  F rdist2;
  double mx,sdx;

  R_rndGen_t(double p,double q, const char* fname)
    : mx(p), sdx(q)
  {
    if ( !strcmp(fname, "rbeta" )) {
        rdist2=rbeta;
    } else if(!strcmp(fname, "rlnorm")) {
        rdist2=rlnorm;
    } else if(!strcmp(fname, "rgamma")) {
        rdist2=rgamma;
    } else if(!strcmp(fname, "runif" )) {
        rdist2=rweibull;
    } else if(!strcmp(fname, "const" )) {
        rdist2=rconst;
    } else {
        error("Undefined random generating function for radii distribution");
    }

  };

  inline double operator()() { return rdist2(mx,sdx); }
};


/**
 * \brief Show aruments of R function call
 *
 * @param args list of arguments
 * @return
 */
SEXP showArgs(SEXP args);

/**
 *
 * @param R_fname
 * @param R_arg
 * @param R_rho
 * @return
 */
SEXP getSingleCall(SEXP R_fname, SEXP R_arg, SEXP R_rho);

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
