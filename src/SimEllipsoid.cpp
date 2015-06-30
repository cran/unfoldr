#define MAX_ITER 100

#include "SimEllipsoid.h"
#include "directions.h"

//static locals
static SEXP spheroid_type_tag;
static int PL = 0;

#define GET_CALL(d,i) getCall( VECTOR_ELT(d->fname,i),VECTOR_ELT(d->args,i),d->rho)
#define GET_NAME(d,i) CHAR(STRING_ELT( VECTOR_ELT(d->fname,i), 0))

SEXP convert_R_EllipsoidSystem( STGM::Spheroids &spheroids);
SEXP convert_R_Ellipses(STGM::IntersectorSpheroidPlaneVec &objects);
STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids);

void _sptr_finalizer(SEXP ext)
{
    checkPtr(ext, spheroid_type_tag);
    STGM::CEllipsoidSystem *sptr = static_cast<STGM::CEllipsoidSystem *>(R_ExternalPtrAddr(ext));
    delete sptr;
    R_ClearExternalPtr(ext);
}

SEXP FinalizeSpheroidSystem(SEXP ext) {
  _sptr_finalizer(ext);
  return R_NilValue;
}

R_Calldata buildRCallSpheroids(SEXP R_param, SEXP R_cond) {
  int nprotect=0;
  R_Calldata d = Calloc(1,R_Calldata_s);

  PROTECT(d->fname = getListElement( R_cond, "rdist"));  ++nprotect;
  if(TYPEOF(d->fname)!=VECSXP) {
      PROTECT(d->args  = getListElement( R_param,"rmulti")); ++nprotect;
      PROTECT(d->rho   = getListElement( R_cond, "rho"  ));  ++nprotect;
      PROTECT(d->call  = getCall(d->fname,d->args,d->rho));  ++nprotect;
  } else {
    PROTECT(d->call = allocVector(VECSXP,3)); ++nprotect;

    /* fname, args are lists of [size, orientation, shape] */
    PROTECT(d->rho   = getListElement( R_cond, "rho"  )); ++nprotect;
    PROTECT(d->args = allocVector(VECSXP,3)); ++nprotect;
    SET_VECTOR_ELT(d->args,0, getListElement( R_param,"size"));
    SET_VECTOR_ELT(d->args,1, getListElement( R_param,"shape"));
    SET_VECTOR_ELT(d->args,2, getListElement( R_param,"orientation"));

    /* size distributions */
    SET_VECTOR_ELT(d->call,0,R_NilValue);
    const char *ftype_size = GET_NAME(d,0);
    if ( !strcmp( ftype_size, "rlnorm") ||
         !strcmp( ftype_size, "rbeta" ) ||
         !strcmp( ftype_size, "rgamma") ||
         !strcmp( ftype_size, "runif" ) ||
         !strcmp( ftype_size, "const" ))
    {;
    } else {
      SET_VECTOR_ELT(d->call,0,GET_CALL(d,0));
    }

    /* shape */
    SET_VECTOR_ELT(d->call,1,R_NilValue);

    /* orientation distributions */
    SET_VECTOR_ELT(d->call,2,R_NilValue);
    const char *ftype_dir = GET_NAME(d,2);
    if ( !strcmp( ftype_dir, "runifdir") ||
         !strcmp( ftype_dir, "rbetaiso" ) ||
         !strcmp( ftype_dir, "rvMisesFisher"))
    {;
    } else {
       SET_VECTOR_ELT(d->call,2,GET_CALL(d,2));
    }

  }

  d->nprotect = nprotect;
  return d;
}

SEXP InitSpheroidSystem(SEXP R_param, SEXP R_cond) {
  int nProtected=0;

  /* get the box */
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));  ++nProtected;
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  SEXP R_spheroidType=R_NilValue;
  PROTECT( R_spheroidType = getListElement( R_cond, "stype" ) );++nProtected;
  const char* stype_str = CHAR( STRING_ELT( R_spheroidType, 0 ));

  STGM::CSpheroid::spheroid_type stype = STGM::CSpheroid::PROLATE;
  if( !strcmp("oblate",stype_str))
    stype = STGM::CSpheroid::OBLATE;

  // print level
  PL = asInteger(getListElement( R_cond,"pl"));

  double lam = asReal(AS_NUMERIC(getListElement( R_param, "lam")));

  /* set up spheroid system */
  STGM::CBox3 box(boxX,boxY,boxZ);

  STGM::CVector3d zaxis(0,0,1);
  STGM::CEllipsoidSystem *sp = new STGM::CEllipsoidSystem(box,lam,zaxis,stype);

  if (!sp) {
      UNPROTECT(nProtected);
      error("Allocation error of spheroid object!\n");
      return R_NilValue;
  }

  SEXP extp;
  spheroid_type_tag = install("SpheroidSystem_TAG");
  PROTECT(extp=R_MakeExternalPtr((void*) sp, spheroid_type_tag, R_NilValue)); ++nProtected;
  R_RegisterCFinalizerEx(extp, (R_CFinalizer_t) _sptr_finalizer, TRUE);

  UNPROTECT(nProtected);
  return extp;
}

SEXP SetupSpheroidSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond)
{
   SEXP R_var, R_ptr;
   PROTECT(R_var = getVar(R_vname,R_env));
   PROTECT(R_ptr = getAttrib(R_var,install("eptr")));

   if(isNull(R_ptr) || isNullPtr(R_ptr,spheroid_type_tag)) {
       R_ptr = InitSpheroidSystem(R_param,R_cond);
       if(PL>100)
         Rprintf("setting pointer to %p \n",R_ExternalPtrAddr(R_ptr));
   }

   STGM::CEllipsoidSystem *sptr = static_cast<STGM::CEllipsoidSystem *>(getExternalPtr(R_ptr));
   sptr->refObjects() = convert_C_Spheroids(R_var);
   setAttrib(R_var, install("eptr"), R_ptr);
   UNPROTECT(2);
   return R_NilValue;
}


SEXP GetEllipsoidSystem(SEXP ext) {
  checkPtr(ext, spheroid_type_tag);

  STGM::CEllipsoidSystem *sp = static_cast<STGM::CEllipsoidSystem *>(getExternalPtr(ext));

  SEXP R_spheroids = R_NilValue;
  STGM::Spheroids &spheroids = sp->refObjects();
  PROTECT(R_spheroids = convert_R_EllipsoidSystem(spheroids));
  setAttrib(R_spheroids, install("eptr"), ext);

  const char *stype=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_spheroids,stype);

  UNPROTECT(1);
  return R_spheroids;

}

SEXP EllipsoidSystem(SEXP R_param, SEXP R_cond) {
  /** Init */
  SEXP ext = PROTECT(InitSpheroidSystem(R_param,R_cond));
  STGM::CEllipsoidSystem *sp = static_cast<STGM::CEllipsoidSystem *>(getExternalPtr(ext));

  R_Calldata call_data = buildRCallSpheroids(R_param,R_cond);

  if(TYPEOF(call_data->fname)!=VECSXP)
    sp->simEllipsoidSysJoint(call_data);
  else
    sp->simEllipsoidSys(call_data);
  deleteRCall(call_data);

  if(PL>100) {
    Rprintf("Simulated %d spheroids.\n",sp->refObjects().size());
  }

  SEXP R_spheroids = R_NilValue;
  if(PL>100) {
      STGM::Spheroids &spheroids = sp->refObjects();
      Rprintf("Convert... \n");
      PROTECT(R_spheroids = convert_R_EllipsoidSystem(spheroids));
  } else {
      PROTECT(R_spheroids = allocVector(VECSXP,0));  /* return empty list */
  }
  setAttrib(R_spheroids, install("eptr"), ext);
  classgets(R_spheroids, getListElement( R_cond, "stype" ));

  UNPROTECT(2);
  return R_spheroids;
}

SEXP IntersectSpheroidSystem(SEXP ext, SEXP R_n, SEXP R_z)
{
  checkPtr(ext, spheroid_type_tag);

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane( n , asReal(R_z));

  STGM::CEllipsoidSystem *sp = static_cast<STGM::CEllipsoidSystem *>(getExternalPtr(ext));
  if(PL>100) Rprintf("Intersect with plane: %d \n", sp->refObjects().size());

  STGM::IntersectorSpheroidPlaneVec objects;
  sp->IntersectWithPlane(objects,plane);

  SEXP R_ellipses;
  PROTECT(R_ellipses = convert_R_Ellipses(objects));

  const char *stype=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_ellipses,stype);

  UNPROTECT(1);
  return R_ellipses;
}

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  SEXP R_tmp, R_ctr, R_u, R_ab, R_angles;

  STGM::Spheroids spheroids;
  for(int i=0; i<length(R_spheroids); i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheroids,i));

      int id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_tmp, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_tmp, "u")));
      PROTECT( R_ab     = AS_NUMERIC( getListElement( R_tmp, "ab")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_tmp, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
      STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

      spheroids.push_back(STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],u,1.0,REAL(R_angles)[0],REAL(R_angles)[1],id));
      UNPROTECT(5);
  }

  return spheroids;
}

/**
 * @brief Convert ellipses to R objects
 *
 * @param objects Intersection object
 * @return R ellipses
 */
SEXP convert_R_Ellipses(STGM::IntersectorSpheroidPlaneVec &objects) {
  int nProtected=0, ncomps=5, dim=2;

  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) ); ++nProtected;

  SEXP names,R_tmp, R_center, R_ab, R_phi, R_id, R_A;

  for(size_t k=0; k<objects.size(); ++k)
  {
      STGM::CEllipse2 &ellipse = objects[k].getEllipse();

      PROTECT(R_tmp = allocVector(VECSXP,ncomps));
      PROTECT(R_id = allocVector(REALSXP, 1));
      PROTECT(R_center = allocVector(REALSXP, dim));
      PROTECT(R_ab = allocVector(REALSXP, dim));
      PROTECT(R_phi = allocVector(REALSXP, 1));
      PROTECT(R_A = allocMatrix(REALSXP,2,2));

      REAL(R_id)[0]=ellipse.Id();
      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];
      REAL(R_ab)[0] = ellipse.a();
      REAL(R_ab)[1] = ellipse.b();
      REAL(R_phi)[0] = ellipse.phi();

      //std::cout << "A: " << ellipse.A << std::endl;

      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          REAL(R_A)[i + dim *j] = ellipse.MatrixA()[i][j];

      SET_VECTOR_ELT(R_tmp,0,R_id);
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_A);
      SET_VECTOR_ELT(R_tmp,3,R_ab);
      SET_VECTOR_ELT(R_tmp,4,R_phi);


      PROTECT(names = allocVector(STRSXP, ncomps));
      SET_STRING_ELT(names, 0, mkChar("id"));
      SET_STRING_ELT(names, 1, mkChar("center"));
      SET_STRING_ELT(names, 2, mkChar("A"));
      SET_STRING_ELT(names, 3, mkChar("ab"));
      SET_STRING_ELT(names, 4, mkChar("phi"));
      setAttrib(R_tmp, R_NamesSymbol, names);

      SET_VECTOR_ELT(R_resultlist,k,R_tmp);
      UNPROTECT(7);
  }

  UNPROTECT(nProtected);
  return(R_resultlist);
}


SEXP convert_R_EllipsoidSystem( STGM::Spheroids &spheroids) {
  int nProtected=0, ncomps=6, dim=3;

  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, spheroids.size()) ); ++nProtected;

  SEXP names, R_center, R_u, R_ab, R_tmp, R_angles, R_rotM, R_id;
  for(size_t k=0;k<spheroids.size();k++)
  {
    STGM::CSpheroid &spheroid = spheroids[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_id = allocVector(REALSXP, 1));
    PROTECT(R_center = allocVector(REALSXP, dim));
    PROTECT(R_u = allocVector(REALSXP, dim));
    PROTECT(R_ab = allocVector(REALSXP,2));
    PROTECT(R_angles = allocVector(REALSXP,2));
    PROTECT(R_rotM = allocMatrix(REALSXP,3,3));

    REAL(R_id)[0]=spheroid.Id();

    STGM::CVector3d m_center(spheroid.Center());
    REAL(R_center)[0]=m_center[0];
    REAL(R_center)[1]=m_center[1];
    REAL(R_center)[2]=m_center[2];

    STGM::CVector3d m_u = spheroid.u();
    REAL(R_u)[0]=m_u[0];
    REAL(R_u)[1]=m_u[1];
    REAL(R_u)[2]=m_u[2];

    REAL(R_ab)[0]=spheroid.a();
    REAL(R_ab)[1]=spheroid.b();

    REAL(R_angles)[0]=spheroid.theta();
    REAL(R_angles)[1]=spheroid.phi();

    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
        REAL(R_rotM)[i + dim *j] = (spheroid.rotationMatrix())[i][j];

    SET_VECTOR_ELT(R_tmp,0,R_id);
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,R_u);
    SET_VECTOR_ELT(R_tmp,3,R_ab);
    SET_VECTOR_ELT(R_tmp,4,R_angles);
    SET_VECTOR_ELT(R_tmp,5,R_rotM);

    PROTECT(names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(names, 0, mkChar("id"));
    SET_STRING_ELT(names, 1, mkChar("center"));
    SET_STRING_ELT(names, 2, mkChar("u"));
    SET_STRING_ELT(names, 3, mkChar("ab"));
    SET_STRING_ELT(names, 4, mkChar("angles"));
    SET_STRING_ELT(names, 5, mkChar("rotM"));
    setAttrib(R_tmp, R_NamesSymbol, names);

    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(8);
  }

  UNPROTECT(nProtected);
  return(R_resultlist);
}

void STGM::CEllipsoidSystem::simEllipsoidSysJoint(R_Calldata d) {
     GetRNGstate();

     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
       num = rpois(m_box.volume()*m_lam);
       ++nTry;
     }
     m_spheroids.reserve(num);

     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double a=0,b=0,shape=1,theta=0, phi=0;
     double *v;

     SEXP Reval=R_NilValue;
     CVector3d u;

     int err = 0;
     for (size_t niter=0; niter<num; niter++)
     {
         //Reval=eval(d->call,d->rho);
         //shape=a/b;
         Reval = R_tryEval(d->call,d->rho,&err);
         if(!err) {
            a=asReal(getListElement(Reval,"a"));
            b=asReal(getListElement(Reval,"b"));
            shape=asReal(getListElement(Reval,"shape"));
            theta=asReal(getListElement(Reval,"theta"));
            phi=asReal(getListElement(Reval,"phi"));

            v=REAL(getListElement(Reval,"u"));
            u[0]=v[0];
            u[1]=v[1];
            u[2]=v[2];

            if(m_stype==CSpheroid::OBLATE)
              std::swap(a,b);

            STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
            m_spheroids.push_back( STGM::CSpheroid(center,a,b,u,shape,theta,phi,m_spheroids.size()+1) );
         } else
           error(_("simEllipsoidSysJoint(): R try error in user defined distribution function."));
     }
     PutRNGstate();
}

void STGM::CEllipsoidSystem::simEllipsoidSys(R_Calldata d)
{
     GetRNGstate();

     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
        num = rpois(m_box.volume()*m_lam);
        ++nTry;
     }
     m_spheroids.reserve(num);

     rdist2_t rdist;
     const char *fname = GET_NAME(d,0);
     if ( !strcmp(fname, "rbeta" )) {
          rdist=&rbeta;
      } else if(!strcmp(fname, "rlnorm")) {
          rdist=&rlnorm;
      } else if(!strcmp(fname, "rgamma")) {
          rdist=&rgamma;
      } else if(!strcmp(fname, "runif" )) {
          rdist=&runif;
      } else if(!strcmp(fname, "const" )) {
          rdist=&rconst;
      }

     const char *fname_dir = GET_NAME(d,2);
     STGM::CSpheroid::direction_type dtype = STGM::CSpheroid::UNIFORM_D;
     if(!strcmp( fname_dir, "rbetaiso" )) {
        dtype=STGM::CSpheroid::BETAISOTROP_D;
     } else if(!strcmp( fname_dir, "rvMisesFisher")) {
        dtype=STGM::CSpheroid::MISES_D;
     }

     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double a=0,b=0,kappa=0,theta=0, phi=0;

     Rboolean size_flag = isNull(VECTOR_ELT(d->call,0));
     Rboolean dir_flag  = isNull(VECTOR_ELT(d->call,1));

     SEXP RCall_size = VECTOR_ELT(d->call,0);
     SEXP RCall_dir  = VECTOR_ELT(d->call,2);

     /* shape is currently only const */
     double shape = asReal(VECTOR_ELT(VECTOR_ELT(d->args,1),0));

     double p1=0,p2=0;
     if(size_flag) {
         p1=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),0));
         if(strcmp(fname, "const" ))
           p2=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),1));

     }

     if(dir_flag)
       kappa = asReal(VECTOR_ELT(VECTOR_ELT(d->args,2),0));

     if(PL>100) {
         Rprintf("Run spheroids  simulation... \n");
         if(size_flag)
           Rprintf("\t size distribution: %s with %f %f \n", fname, p1,p2);
         else
           showArgs(VECTOR_ELT(d->args,0)); // size
         Rprintf("\t constant shape parameter: %f \n", shape);
         if(dir_flag)
          Rprintf("\t directional distribution: %s  with %f \n", fname_dir, kappa);
         else
           showArgs(VECTOR_ELT(d->args,2)); // size
     }

     /* loop over all */
     double *v;
     CVector3d u;
     for (size_t niter=0; niter<num; niter++)
     {
         if(size_flag) {
             b = rdist(p1,p2);
         } else {
             b=asReal(eval(RCall_size,d->rho));
         }
         a = b * shape;
         if(m_stype==CSpheroid::OBLATE) std::swap(a,b);

         /* direction */
         if(dir_flag) {
           switch(dtype) {
             case 0:
               runidir(u.ptr(),theta,phi); break;
             case 1:
               if(kappa<1e-8) {
                 u = (runif(0.0,1.0)<0.5) ? m_mu : -m_mu;
               } else {
                 rOhserSchladitz(u.ptr(),m_mu.ptr(), kappa, theta, phi);
               }
               break;
             case 2:
               if(kappa<1e-8)
                 runidir(u.ptr(),theta,phi);
               else
                 rVonMisesFisher(u.ptr(), m_mu.ptr(), kappa, theta, phi);
               break;
           }
         } else {
             /**
              * @todo !!! test code !!!
              * */
             v = NUMERIC_POINTER(eval(RCall_dir,d->rho));
             u[0]=v[0]; u[1]=v[1]; u[2]=v[2];
         }

        STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
        m_spheroids.push_back( STGM::CSpheroid(center,a,b,u,shape,theta,phi,m_spheroids.size()+1) );
     }
     PutRNGstate();
}

void STGM::CEllipsoidSystem::IntersectWithPlane(STGM::IntersectorSpheroidPlaneVec &objects,STGM::CPlane &plane)
{
  /// Intersect only objects fully inside the observation window
    for(size_t i=0; i<m_spheroids.size(); ++i) {
         STGM::IntersectorSpheroid intersector( m_spheroids[i], plane, m_box.m_size);
         if(intersector.FindIntersection()) {
             //if(intersector.getEllipse().isInWindow(m_win))
               objects.push_back( intersector );
         }
     }
}
