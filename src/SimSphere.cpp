/**
 * simSphere.cpp
 *
 *  Created on: 07.05.2015
 *      Author: franke
 */

#include "SimSphere.h"

//static locals
static SEXP sphere_type_tag;
static int PL = 0;

#define MAX_ITER 100

SEXP convert_R_SphereSystem(STGM::Spheres& spheres);
SEXP convert_R_Circles(STGM::IntersectorSpherePlaneVec& objects);
STGM::Spheres convert_C_Spheres(SEXP R_spheres);

R_Calldata buildRCallSpheres(SEXP R_param,SEXP R_cond) {
  int nprotect=0;
  R_Calldata d = Calloc(1,R_Calldata_s);
  d->call = R_NilValue;
  PROTECT(d->fname = getListElement( R_cond, "rdist")); ++nprotect;
  PROTECT(d->rho   = getListElement( R_cond, "rho"  )); ++nprotect;
  PROTECT(d->args  = getListElement( R_param,"radii")); ++nprotect;

  /* radii distribution */
  const char *ftype = CHAR(STRING_ELT(d->fname, 0));
  if ( !strcmp( ftype, "rlnorm") ||
       !strcmp( ftype, "rbeta" ) ||
       !strcmp( ftype, "rgamma") ||
       !strcmp( ftype, "runif" ) ||
       !strcmp( ftype, "const" ))
  {;
  } else {
    PROTECT(d->call = getCall(d->fname,d->args,d->rho)); ++nprotect;
  }
  d->nprotect = nprotect;
  return d;
}

void _sphere_finalizer(SEXP ext)
{
    checkPtr(ext, sphere_type_tag);
    STGM::CBoolSphereSystem *sptr = static_cast<STGM::CBoolSphereSystem *>(R_ExternalPtrAddr(ext));
    delete sptr;
    R_ClearExternalPtr(ext);
}


SEXP FinalizeSphereSystem(SEXP ext) {
  _sphere_finalizer(ext);
  return R_NilValue;
}

SEXP InitSphereSystem(SEXP R_param, SEXP R_cond) {
  int nProtected=0;

  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));  ++nProtected;
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  double lam   = asReal(AS_NUMERIC( getListElement( R_param, "lam")));
  // print level
  PL = asInteger(getListElement( R_cond,"pl"));

  STGM::CBox3 box(boxX,boxY,boxZ);

  /* set up sphere system */
  STGM::CBoolSphereSystem *sp = new STGM::CBoolSphereSystem(box,lam);

  SEXP extp;
  sphere_type_tag = install("SphereSystem_TAG");
  PROTECT(extp = R_MakeExternalPtr((void*) sp, sphere_type_tag, R_NilValue)); ++nProtected;
  R_RegisterCFinalizerEx(extp, (R_CFinalizer_t) _sphere_finalizer, TRUE);

  UNPROTECT(nProtected);
  return extp;
}


template< typename  F>
void STGM::CBoolSphereSystem::simSpheres(F f) {
  double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
         m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
         m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

  /* loop over all */
  for (size_t niter=0;niter<num; niter++) {
      STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
      m_spheres.push_back( STGM::CSphere(center, f(), m_spheres.size()+1));
  }

}

void STGM::CBoolSphereSystem::simSphereSys(R_Calldata d) {
   int nTry = 0;

   GetRNGstate();
   while(num==0 && nTry<MAX_ITER) {
      num = rpois(m_box.volume()*m_lam);
      ++nTry;
   }
   m_spheres.reserve(num);

   if(isNull(d->call)) {
       /* get arguments */
       double p1=REAL_ARG_LIST(d->args,0),p2=0;
       const char *fname = CHAR(STRING_ELT(d->fname,0));
       if(strcmp(fname, "const" ))
         p2=REAL_ARG_LIST(d->args,1); /* ACHTUNG: 'const' function braucht 2 Argumente */

       rdist2_t rdist2=rconst;
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

       R_rndGen_t<rdist2_t> rrandom(p1,p2,rdist2);

       /* simulate with R's random generating functions */
       simSpheres<R_rndGen_t<rdist2_t> >(rrandom);

   } else {
       /* eval R call */
       R_eval_t<double> reval(d->call, d->rho);
       simSpheres<R_eval_t<double> &>(reval);

   }
   PutRNGstate();
}


SEXP GetSphereSystem(SEXP ext)
{
  checkPtr(ext, sphere_type_tag);

  STGM::CBoolSphereSystem *sp = static_cast<STGM::CBoolSphereSystem *>(getExternalPtr(ext));
  SEXP R_spheres = R_NilValue;
  STGM::Spheres &spheres = sp->refObjects();
  PROTECT(R_spheres = convert_R_SphereSystem(spheres));
  setAttrib(R_spheres, install("eptr"), ext);

  SET_CLASS_NAME(R_spheres,"spheres");

  UNPROTECT(1);
  return R_spheres;
}

SEXP SetupSphereSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond)
{
  SEXP R_var, R_ptr;
  PROTECT(R_var = getVar(R_vname,R_env));
  PROTECT(R_ptr = getAttrib(R_var,install("eptr")));

  if(isNull(R_ptr) || isNullPtr(R_ptr,sphere_type_tag)) {
     R_ptr = InitSphereSystem(R_param,R_cond);
     if(PL>100)
       Rprintf("setting pointer to %p \n",R_ExternalPtrAddr(R_ptr));
  }

  STGM::CBoolSphereSystem *sp = static_cast<STGM::CBoolSphereSystem *>(getExternalPtr(R_ptr));
  sp->refObjects() = convert_C_Spheres(R_var);

  setAttrib(R_var, install("eptr"), R_ptr);
  UNPROTECT(2);
  return R_NilValue;
}

SEXP SphereSystem(SEXP R_param, SEXP R_cond) {
  SEXP ext = PROTECT(InitSphereSystem(R_param,R_cond));
  STGM::CBoolSphereSystem *sp = static_cast<STGM::CBoolSphereSystem *>(getExternalPtr(ext));

  R_Calldata cdata = buildRCallSpheres(R_param,R_cond);

  if(PL>100) Rprintf("Simulate... \n");
  sp->simSphereSys(cdata);
  deleteRCall(cdata);

  if(PL>100) Rprintf("Simulated %d spheres.\n", sp->getNumSpheres());

  SEXP R_spheres=R_NilValue;
  if(PL>100) {
    Rprintf("Convert... \n");
    STGM::Spheres &spheres = sp->refObjects();
    PROTECT(R_spheres = convert_R_SphereSystem(spheres));
  } else {
    PROTECT(R_spheres = allocVector(VECSXP,0));  /* return empty list */
  }
  setAttrib(R_spheres, install("eptr"), ext);
  SET_CLASS_NAME(R_spheres,"spheres");

  UNPROTECT(2);
  return R_spheres;
}

SEXP IntersectSphereSystem(SEXP ext, SEXP R_n, SEXP R_z) {
  checkPtr(ext, sphere_type_tag);

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane( n , asReal(R_z));

  STGM::CBoolSphereSystem *sp = static_cast<STGM::CBoolSphereSystem *>(getExternalPtr(ext));
  if(PL>100) Rprintf("Intersect with plane: %d \n", sp->refObjects().size());

  STGM::IntersectorSpherePlaneVec objects;
  sp->IntersectWithPlane(objects,plane);

  return convert_R_Circles(objects);
}


void STGM::CBoolSphereSystem::IntersectWithPlane(STGM::IntersectorSpherePlaneVec &objects, STGM::CPlane &plane)
{
  /// Intersect only objects fully inside the observation window
   for(size_t i=0; i<m_spheres.size(); ++i) {
        STGM::IntersectorSphere intersector( m_spheres[i], plane, m_box.m_size);
        if(intersector.FindIntersection()) {
          //if(intersector.getCircle().isInWindow(m_win))
              objects.push_back( intersector );
        }
    }
}

SEXP convert_R_SphereSystem(STGM::Spheres& spheres) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, spheres.size()) );

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<spheres.size();k++)
  {
    STGM::CSphere &sphere = spheres[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, 3));

    REAL(R_center)[0]=sphere.center()[0];
    REAL(R_center)[1]=sphere.center()[1];
    REAL(R_center)[2]=sphere.center()[2];

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

    PROTECT(R_names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(R_names, 0, mkChar("id"));
    SET_STRING_ELT(R_names, 1, mkChar("center"));
    SET_STRING_ELT(R_names, 2, mkChar("r"));

    setAttrib(R_tmp, R_NamesSymbol, R_names);
    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(3);
  }

  UNPROTECT(1);
  return R_resultlist;
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres) {
  SEXP R_tmp, R_ctr;
  int id=0,
      N=length(R_spheres);
  STGM::Spheres spheres;
  spheres.reserve(N);

  double r;
  for(int i=0; i<N; i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheres,i));
      PROTECT(R_ctr = AS_NUMERIC( getListElement( R_tmp, "center")));
      id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
      r = asReal(getListElement( R_tmp, "r"));

      spheres.push_back(STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id));
      UNPROTECT(2);
  }

  return spheres;
}


SEXP convert_R_Circles(STGM::IntersectorSpherePlaneVec& objects) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<objects.size();k++)
  {
     STGM::CCircle3 &circle = objects[k].getCircle();

     PROTECT(R_tmp = allocVector(VECSXP,ncomps));
     PROTECT(R_center = allocVector(REALSXP, 3));

     REAL(R_center)[0]=circle.center()[0];
     REAL(R_center)[1]=circle.center()[1];
     REAL(R_center)[2]=circle.center()[2];

     SET_VECTOR_ELT(R_tmp,0,ScalarInteger(circle.Id()));
     SET_VECTOR_ELT(R_tmp,1,R_center);
     SET_VECTOR_ELT(R_tmp,2,ScalarReal(circle.r()));

     PROTECT(R_names = allocVector(STRSXP, ncomps));
     SET_STRING_ELT(R_names, 0, mkChar("id"));
     SET_STRING_ELT(R_names, 1, mkChar("center"));
     SET_STRING_ELT(R_names, 2, mkChar("r"));

     setAttrib(R_tmp, R_NamesSymbol, R_names);
     SET_VECTOR_ELT(R_resultlist,k,R_tmp);
     UNPROTECT(3);
   }

   UNPROTECT(1);
   return R_resultlist;
 }
