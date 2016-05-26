/**
 *  @file SimEllipsoid.cpp
 *
 *  @date: 10.01.2014
 *  @author: M. Baaske
 */

#define MAX_ITER 100

#include "unfold.h"
#include "directions.h"
#include "SimEllipsoid.h"


//static locals
static SEXP spheroid_type_tag;
static int PL = 0;

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_j+DIM*_i] = (M)[_i][_j];         \
} while(0)

#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_j+DIM*_i];         \
} while(0)

#define GET_OBJECT_CLASS(RS) translateChar(asChar(getAttrib( (RS), R_ClassSymbol)))

SEXP convert_R_EllipsoidSystem( STGM::Spheroids &spheroids, STGM::CBox3 &box);
SEXP convert_R_Ellipses_all(STGM::IntersectorSpheroids &objects);
SEXP convert_R_Ellipses_trunc(STGM::IntersectorSpheroids &objects);

extern STGM::CSphere convert_C_Sphere(SEXP R_sphere);
extern STGM::CCylinder convert_C_Cylinder(SEXP R_cyl);
extern STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid);

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids);

void _free_spheroids(STGM::CEllipsoidSystem *sp){
  if(NULL==sp) return;
  //Rprintf("free spheroids... \n");
  sp->~CEllipsoidSystem();
  Free(sp);
}

void _sptr_finalizer(SEXP ext)
{
  //if (NULL == R_ExternalPtrAddr(ext))
  //    return;
  //char *sp = (char *) R_ExternalPtrAddr(ext);
  //Rprintf("finalize... \n");

  checkPtr(ext, spheroid_type_tag);
  STGM::CEllipsoidSystem *sp = (STGM::CEllipsoidSystem *)(R_ExternalPtrAddr(ext));
  _free_spheroids(sp);
  R_ClearExternalPtr(ext);
}

SEXP FinalizeSpheroidSystem(SEXP ext) {
  _sptr_finalizer(ext);
  return R_NilValue;
}



SEXP CreateExternalPointer(STGM::CEllipsoidSystem *sp) {
  spheroid_type_tag = install("SpheroidSystem_TAG");
  SEXP Rextp=R_NilValue;
  PROTECT(Rextp=R_MakeExternalPtr((void*) sp, spheroid_type_tag, R_NilValue));
  R_RegisterCFinalizerEx(Rextp, (R_CFinalizer_t) _sptr_finalizer, TRUE);

  UNPROTECT(1);
  return Rextp;
}

STGM::CEllipsoidSystem * InitSpheroidSystem(SEXP R_param, SEXP R_cond) {
  /* get the box */
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  SEXP R_spheroidType = R_NilValue;
  PROTECT( R_spheroidType = getListElement( R_cond, "stype" ) );
  const char* stype_str = CHAR( STRING_ELT( R_spheroidType, 0 ));

  STGM::CSpheroid::spheroid_type stype = STGM::CSpheroid::PROLATE;
  if( !strcmp("oblate",stype_str))
    stype = STGM::CSpheroid::OBLATE;

  // set print level
  PL = asInteger(getListElement( R_cond,"pl"));
  double lam = asReal(AS_NUMERIC(getListElement( R_param, "lam")));

  /* set up spheroid system */
  STGM::CBox3 box(boxX,boxY,boxZ);

  double *mu = NUMERIC_POINTER( getListElement( R_cond, "mu"));
  STGM::CVector3d maxis(mu[0],mu[1],mu[2]);

  STGM::CEllipsoidSystem *sp = (STGM::CEllipsoidSystem *)Calloc(1,STGM::CEllipsoidSystem);

  try {
      new(sp)STGM::CEllipsoidSystem(box,lam,maxis,stype);
  } catch(...) {
      error(_("InitSpheroidSystem(): Allocation error."));
  }

  UNPROTECT(2);
  return sp;
}

SEXP SetupSpheroidSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond)
{
   SEXP R_var, R_ptr;
   PROTECT(R_var = getVar(R_vname,R_env));
   PROTECT(R_ptr = getAttrib(R_var,install("eptr")));

   if(isNull(R_ptr) || isNullPtr(R_ptr,spheroid_type_tag)) {
       R_ptr = CreateExternalPointer( InitSpheroidSystem(R_param,R_cond) );
       if(PL>100) Rprintf("setting pointer to %p \n",R_ExternalPtrAddr(R_ptr));
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
  PROTECT(R_spheroids = convert_R_EllipsoidSystem(spheroids,sp->box()));
  setAttrib(R_spheroids, install("eptr"), ext);

  const char *stype=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_spheroids,stype);

  UNPROTECT(1);
  return R_spheroids;

}

SEXP EllipsoidSystem(SEXP R_param, SEXP R_cond) {
  /** Init */
  STGM::CEllipsoidSystem *sp = InitSpheroidSystem(R_param,R_cond);
  R_Calldata call_data = getRCallParam(R_param,R_cond);

  if(TYPEOF(call_data->fname) != VECSXP) {
      sp->simSysJoint(call_data);
  } else {
      //const char *ftype_size = GET_NAME(call_data,0);
      if(!strcmp(GET_NAME(call_data,0), "rbinorm") ) {
          sp->simBivariate(call_data);
      } else {
          sp->simEllipsoidSys(call_data);
      }
  }
  deleteRCall(call_data);

  if(PL>100) {
    Rprintf("Simulated %d spheroids: %p \n",sp->refObjects().size(), sp);
  }

  SEXP R_spheroids = R_NilValue;
  if(PL>100) {
      STGM::Spheroids &spheroids = sp->refObjects();
      Rprintf("Convert... \n");
      PROTECT(R_spheroids = convert_R_EllipsoidSystem(spheroids,sp->box()));
  } else {
      PROTECT(R_spheroids = allocVector(VECSXP,0));  /* return empty list */
  }

  SEXP Rextp=R_NilValue;
  PROTECT(Rextp=CreateExternalPointer(sp));

  setAttrib(R_spheroids, install("eptr"), Rextp);
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
  if(PL>100) {
    Rprintf("Intersect with plane: %d , %p \n", sp->refObjects().size(), sp);
  }

  STGM::IntersectorSpheroids objects;
  sp->IntersectWithPlane(objects,plane);

  SEXP R_ellipses;
  if(PL==10) {
    PROTECT(R_ellipses = convert_R_Ellipses_trunc(objects));
  } else {
    PROTECT(R_ellipses = convert_R_Ellipses_all(objects));
  }

  const char *stype=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_ellipses,stype);

  UNPROTECT(1);
  return R_ellipses;
}


SEXP SimulateSpheroidsAndIntersect(SEXP R_param, SEXP R_cond) {
  /** Init */
  STGM::CEllipsoidSystem *sp = InitSpheroidSystem(R_param,R_cond);
  R_Calldata call_data = getRCallParam(R_param,R_cond);

  if(TYPEOF(call_data->fname)!=VECSXP) {
    sp->simSysJoint(call_data);
  } else {
      //const char *ftype_size = GET_NAME(call_data,0);
      if(!strcmp(GET_NAME(call_data,0), "rbinorm") ) {
        sp->simBivariate(call_data);
      } else {
        sp->simEllipsoidSys(call_data);
      }
  }

  deleteRCall(call_data);

  double dz = asReal(getListElement(R_cond,"dz"));
  STGM::CPlane plane(STGM::CVector3d(1,0,0),dz);

  STGM::IntersectorSpheroids objects;
  sp->IntersectWithPlane(objects,plane);

  SEXP R_ellipses;
  if(PL==10) {
    PROTECT(R_ellipses = convert_R_Ellipses_trunc(objects));
  } else {
    PROTECT(R_ellipses = convert_R_Ellipses_all(objects));
  }
  const char *stype=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_ellipses,stype);

  _free_spheroids(sp);

  UNPROTECT(1);
  return R_ellipses;

}

SEXP convert_C2R_ellipses(STGM::Ellipses &ellipses) {
  int nProtected=0, dim=2, ncomps=8;
  int n = ellipses.size();

  SEXP names, R_resultlist;
  PROTECT(names = allocVector(STRSXP, ncomps));   ++nProtected;
  PROTECT(R_resultlist = allocVector(VECSXP,n));  ++nProtected;

  SET_STRING_ELT(names, 0, mkChar("id"));
  SET_STRING_ELT(names, 1, mkChar("center"));
  SET_STRING_ELT(names, 2, mkChar("ab"));
  SET_STRING_ELT(names, 3, mkChar("minor"));
  SET_STRING_ELT(names, 4, mkChar("major"));
  SET_STRING_ELT(names, 5, mkChar("A"));
  SET_STRING_ELT(names, 6, mkChar("phi"));
  SET_STRING_ELT(names, 7, mkChar("rot"));

  SEXP R_tmp,R_minor,R_major,R_A,R_center,R_ab;

  for(int i = 0; i < n; i++) {
      STGM::CEllipse2 & ellipse = ellipses[i];
      PROTECT(R_tmp = allocVector(VECSXP,ncomps));
      PROTECT(R_center = allocVector(REALSXP, dim));
      PROTECT(R_ab = allocVector(REALSXP, dim));
      PROTECT(R_minor = allocVector(REALSXP, dim));
      PROTECT(R_major = allocVector(REALSXP, dim));
      PROTECT(R_A = allocMatrix(REALSXP, dim,dim));

      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];

      REAL(R_minor)[0] = ellipse.minorAxis()[0];
      REAL(R_minor)[1] = ellipse.minorAxis()[1];

      REAL(R_major)[0] = ellipse.majorAxis()[0];
      REAL(R_major)[1] = ellipse.majorAxis()[1];

      REAL(R_ab)[0] = ellipse.a();
      REAL(R_ab)[1] = ellipse.b();

      COPY_C2R_MATRIX(ellipse.MatrixA(),R_A,dim);

      setAttrib(R_tmp, R_NamesSymbol, names);
      SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_ab);
      SET_VECTOR_ELT(R_tmp,3,R_minor);
      SET_VECTOR_ELT(R_tmp,4,R_major);
      SET_VECTOR_ELT(R_tmp,5,R_A);
      SET_VECTOR_ELT(R_tmp,6,ScalarReal(ellipse.phi()));
      SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.rot()));

      SET_VECTOR_ELT(R_resultlist,i,R_tmp);
      UNPROTECT(6);
  }

  UNPROTECT(nProtected);
  return R_resultlist;
}

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  SEXP R_tmp, R_ctr, R_u, R_ab, R_angles;

  SEXP R_label = R_NilValue;
  PROTECT(R_label = getAttrib(R_spheroids, install("label")));
  const char *label = "N";
  if(!isNull(R_label)) {
    label = translateChar(asChar(R_label));
  }

  int interior = 1;
  double radius = 0;
  STGM::Spheroids spheroids;
  for(int i=0; i<length(R_spheroids); i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheroids,i));

      if(!isNull(getAttrib(R_tmp, install("label"))))
        label = translateChar(asChar(getAttrib(R_tmp, install("label"))));

      if(!isNull(getAttrib(R_tmp, install("interior"))))
        interior = asLogical(getAttrib(R_tmp, install("interior")));

      if(!isNull(getAttrib(R_tmp, install("radius"))))
        radius = asReal(getAttrib(R_tmp, install("radius")));

      int id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_tmp, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_tmp, "u")));
      PROTECT( R_ab     = AS_NUMERIC( getListElement( R_tmp, "ab")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_tmp, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
      STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

      spheroids.push_back(STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],u,1.0,REAL(R_angles)[0],REAL(R_angles)[1],radius,id,label,interior));
      UNPROTECT(5);

  }

  UNPROTECT(1);
  return spheroids;
}


SEXP convert_R_Ellipses_trunc(STGM::IntersectorSpheroids &objects) {
  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp=R_NilValue;
  for(size_t k=0; k<objects.size(); ++k)   {
      STGM::CEllipse2 &ellipse = objects[k].getEllipse();
      PROTECT(R_tmp = allocVector(VECSXP,4));
      SET_VECTOR_ELT(R_tmp,0,ScalarReal(ellipse.a()));                   /* size C */
      SET_VECTOR_ELT(R_tmp,1,ScalarReal(ellipse.b()));                   /* size A */
      SET_VECTOR_ELT(R_tmp,2,ScalarReal(ellipse.a()/ellipse.b()));       /* shape  */
      SET_VECTOR_ELT(R_tmp,3,ScalarReal(ellipse.phi()));                 /* orientation, here in [0,2pi] */

      SET_VECTOR_ELT(R_resultlist,k,R_tmp);
      UNPROTECT(1);
  }
  UNPROTECT(1);
  return(R_resultlist);
}

/**
 * @brief Convert ellipses to R objects
 *
 * @param objects Intersection object
 * @return R ellipses
 */
SEXP convert_R_Ellipses_all(STGM::IntersectorSpheroids &objects) {
  int nProtected=0, ncomps=5, dim=2;

  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) ); ++nProtected;

  SEXP names,R_tmp,R_center,R_ab,R_A;

  for(size_t k=0; k<objects.size(); ++k)
  {
      STGM::CEllipse2 &ellipse = objects[k].getEllipse();

      PROTECT(R_tmp = allocVector(VECSXP,ncomps));
      PROTECT(R_center = allocVector(REALSXP, dim));
      PROTECT(R_ab = allocVector(REALSXP, dim));
      PROTECT(R_A = allocMatrix(REALSXP,2,2));

      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];
      REAL(R_ab)[0] = ellipse.a();
      REAL(R_ab)[1] = ellipse.b();

      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          REAL(R_A)[i + dim *j] = ellipse.MatrixA()[i][j];

      SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_A);
      SET_VECTOR_ELT(R_tmp,3,R_ab);
      SET_VECTOR_ELT(R_tmp,4,ScalarReal(ellipse.phi()));

      PROTECT(names = allocVector(STRSXP, ncomps));
      SET_STRING_ELT(names, 0, mkChar("id"));
      SET_STRING_ELT(names, 1, mkChar("center"));
      SET_STRING_ELT(names, 2, mkChar("A"));
      SET_STRING_ELT(names, 3, mkChar("ab"));
      SET_STRING_ELT(names, 4, mkChar("phi"));
      setAttrib(R_tmp, R_NamesSymbol, names);

      SET_VECTOR_ELT(R_resultlist,k,R_tmp);
      UNPROTECT(5);
  }

  UNPROTECT(nProtected);
  return(R_resultlist);
}

STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  SEXP R_ctr, R_u, R_ab, R_angles;

  int id = asInteger (AS_INTEGER( getListElement( R_spheroid, "id")));
  PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_spheroid, "center")));
  PROTECT( R_u      = AS_NUMERIC( getListElement( R_spheroid, "u")));
  PROTECT( R_ab     = AS_NUMERIC( getListElement( R_spheroid, "ab")));
  PROTECT( R_angles = AS_NUMERIC( getListElement( R_spheroid, "angles")));

  STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
  STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

  int interior = 1;
  double radius = 0;
  const char *label = "N";
  if(!isNull(getAttrib(R_spheroid, install("label"))))
    label = translateChar(asChar(getAttrib(R_spheroid, install("label"))));

  if(!isNull(getAttrib(R_spheroid, install("interior"))))
    interior = asLogical(getAttrib(R_spheroid, install("interior")));
  if(!isNull(getAttrib(R_spheroid, install("radius"))))
    radius = asReal(getAttrib(R_spheroid, install("radius")));

  UNPROTECT(4);
  return STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],u,1.0,REAL(R_angles)[0],REAL(R_angles)[1],radius,id,label,interior);
}


STGM::Ellipses convert_C_Ellipses(SEXP R_ellipses)
{
  int id=0;
  STGM::Ellipses ellipses;
  SEXP R_tmp, R_ctr,R_ab, R_minor, R_major, R_A;

  double rot = 0;
  for(int i=0; i<length(R_ellipses); i++) {
     PROTECT(R_tmp   = VECTOR_ELT(R_ellipses,i));
     PROTECT(R_ctr   = AS_NUMERIC( getListElement( R_tmp, "center")));
     PROTECT(R_A     = AS_NUMERIC( getListElement( R_tmp, "A")));
     PROTECT(R_minor = AS_NUMERIC( getListElement( R_tmp, "minor")));
     PROTECT(R_major = AS_NUMERIC( getListElement( R_tmp, "major")));
     PROTECT(R_ab    = AS_NUMERIC( getListElement( R_tmp, "ab")));

     id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
     rot= asReal(AS_NUMERIC(getListElement( R_tmp, "rot")));

     /**
      * BUG: Constructor with matrix A has a bug to determine the correct angle phi
      */

     STGM::CPoint2d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1]);
     STGM::CPoint2d minorAxis(REAL(R_minor)[0],REAL(R_minor)[1]);
     STGM::CPoint2d majorAxis(REAL(R_major)[0],REAL(R_major)[1]);

     ellipses.push_back(STGM::CEllipse2(ctr,majorAxis,minorAxis,REAL(R_ab)[0],REAL(R_ab)[1],id,rot));
     UNPROTECT(6);
  }
  return ellipses;
}


SEXP UpdateIntersections(SEXP Rs, SEXP R_box) {
    int nProtected=0;

    SEXP R_BoxX, R_BoxY, R_BoxZ;
    PROTECT( R_BoxX = AS_NUMERIC( getListElement( R_box, "xrange" ) ) ); ++nProtected;
    PROTECT( R_BoxY = AS_NUMERIC( getListElement( R_box, "yrange" ) ) ); ++nProtected;
    PROTECT( R_BoxZ = AS_NUMERIC( getListElement( R_box, "zrange" ) ) ); ++nProtected;

    STGM::CBox3 box(REAL(R_BoxX)[1],REAL(R_BoxY)[1],REAL(R_BoxZ)[1]);
    const std::vector<STGM::CPlane> &planes = box.getPlanes();

    SEXP R_ret = R_NilValue;
    PROTECT(R_ret = allocVector(INTSXP,length(Rs)));
    ++nProtected;

    const char * name = GET_OBJECT_CLASS(Rs);
    if( !strcmp(name, "prolate" ) || !strcmp(name, "oblate" ) || !strcmp(name, "spheroid" )) {
        for(int k=0;k<length(Rs);k++) {
                STGM::CSpheroid sp = convert_C_Spheroid(VECTOR_ELT(Rs,k));
                STGM::IntersectorSpheroid intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!strcmp(name, "cylinder" )) {
        for(int k=0;k<length(Rs);k++) {
                STGM::CCylinder sp = convert_C_Cylinder(VECTOR_ELT(Rs,k));
                STGM::IntersectorCylinder intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!strcmp(name, "sphere" )) {
      for(int k=0;k<length(Rs);k++) {
                STGM::CSphere sp = convert_C_Sphere(VECTOR_ELT(Rs,k));
                STGM::IntersectorSphere intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
      }
    } else {
        error("Unknown class object.");
    }

    UNPROTECT(nProtected);
    return R_ret;

  }


SEXP GetMaxRadius(SEXP ext)
{
  checkPtr(ext, spheroid_type_tag);
  STGM::CEllipsoidSystem *sp = static_cast<STGM::CEllipsoidSystem *>(getExternalPtr(ext));
  return ScalarReal(sp->maxR());
}

SEXP convert_R_EllipsoidSystem( STGM::Spheroids &spheroids, STGM::CBox3 &box) {
  int nProtected=0, ncomps=6, dim=3;

  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, spheroids.size()) ); ++nProtected;

  SEXP names, R_center, R_u, R_ab, R_tmp, R_angles, R_rotM;
  // get lateral bounding planes
  const STGM::LateralPlanes &planes = box.getLateralPlanes();

  for(size_t k=0;k<spheroids.size();k++)
  {
    STGM::CSpheroid &spheroid = spheroids[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, dim));
    PROTECT(R_u = allocVector(REALSXP, dim));
    PROTECT(R_ab = allocVector(REALSXP,2));
    PROTECT(R_angles = allocVector(REALSXP,2));
    PROTECT(R_rotM = allocMatrix(REALSXP,3,3));

    // projection
    STGM::CEllipse2 ellipse = spheroid.spheroidProjection();
    // check intersection
    STGM::IntersectorSpheroid intersector(spheroid , box.m_size );
    Rboolean interior = (Rboolean) TRUE;
    for(size_t j=0; j<planes.size() ; ++j) {
         if( intersector(planes[j])) {
           interior = (Rboolean) FALSE;
           break;
         }
    }

    const STGM::CVector3d &m_center = spheroid.center();
    REAL(R_center)[0]=m_center[0];
    REAL(R_center)[1]=m_center[1];
    REAL(R_center)[2]=m_center[2];

    const STGM::CVector3d &m_u = spheroid.u();
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

    //STGM::CMatrix3d &M = spheroid.rotationMatrix();
    //COPY_C2R_MATRIX(M,R_rotM,dim);

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(spheroid.Id()));
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
    setAttrib(R_tmp, install("radius"), ScalarReal(spheroid.radius()));
    setAttrib(R_tmp, install("label"), mkString(spheroid.label()) );
    setAttrib(R_tmp, install("interior"), ScalarLogical(interior));
    setAttrib(R_tmp, install("area"), ScalarReal(ellipse.area()));

    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(7);
  }

  UNPROTECT(nProtected);
  return(R_resultlist);
}

void STGM::CEllipsoidSystem::simSysJoint(R_Calldata d) {
     GetRNGstate();

     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
       num = rpois(m_box.volume()*m_lam);
       ++nTry;
     }
     m_spheroids.reserve(num);
     // set spheroid label
     const char *label = translateChar(asChar(d->label));
     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double a=0,b=0,shape=1,theta=0, phi=0,r=0;
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
            m_spheroids.push_back( STGM::CSpheroid(center,a,b,u,shape,theta,phi,r,m_spheroids.size()+1,label) );
         } else
           error(_("simSysJoint(): R try error in user defined distribution function."));
     }
     PutRNGstate();
}

/**
 * @brief Bivariate size-shape distribution,
 *         the major semi-axis is logN distributed
 *
 * @param d data of R call (no R call used here)
 */
void STGM::CEllipsoidSystem::simBivariate(R_Calldata d) {
      GetRNGstate();

      double mx=asReal(getListElement(VECTOR_ELT( d->args, 0),"mx"));
      double my=asReal(getListElement(VECTOR_ELT( d->args, 0),"my"));
      double sdx=asReal(getListElement(VECTOR_ELT( d->args, 0),"sdx"));
      double sdy=asReal(getListElement(VECTOR_ELT( d->args, 0),"sdy"));
      double rho=asReal(getListElement(VECTOR_ELT( d->args, 0),"rho"));
      double kappa = asReal(getListElement(VECTOR_ELT(d->args,2),"kappa"));

      /* get Poisson parameter */
      double p[4], sdx2 = SQR(sdx), mu=0;
      // set spheroid label
      const char *label = translateChar(asChar(d->label));
      // cumulative probabilities
      cum_prob_k(mx,sdx2,m_box.m_up[0],m_box.m_up[1],m_box.m_up[2],p,&mu);

      if(PL>100) {
         Rprintf("Spheroids (perfect) simulation, bivariate lognormal length/shape: \n");
         Rprintf("\t size distribution: %s with %f %f %f %f %f\n", GET_NAME(d,0), mx,my,sdx,sdy,rho);
         Rprintf("\t directional distribution: %s  with %f \n", GET_NAME(d,2), kappa);
         Rprintf("\t cum sum of probabilities: %f, %f, %f, %f \n",p[0],p[1],p[2],p[3]);
         Rprintf("\t set label: %s to character: \n",label);
      }
      int nTry=0;
      while(num==0 && nTry<MAX_ITER) {
         num = rpois(mu*m_lam);
         ++nTry;
      }
      m_spheroids.reserve(num);

      CVector3d u;
      double x,y,a,s,
             phi,theta;

      int k=0;
      double r=0;
      for (size_t niter=0; niter<num; niter++)
      {
          /* sample major semi-axis a, shorter semi-axis is: c=a*s */
          rbinorm(mx,sdx,my,sdy,rho,x,y);
          s=1.0/(1.0+exp(-y));
          a=exp(x);

          /* sample orientation */
          if(kappa<1e-8)
            u = (runif(0.0,1.0)<0.5) ? m_mu : -m_mu;
          else rOhserSchladitz(u.ptr(),m_mu.ptr(),kappa,theta,phi);

          /* sample positions conditionally of radii distribution */
          sample_k(p,&k);
          r=rlnorm(mx+k*sdx2,sdx);
          if(m_maxR<r) m_maxR=r;

          if(!R_FINITE(r))
            warning(_("simEllipsoidSysBivariat(): Some NA/NaN, +/-Inf produced"));

          STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
                                 runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
                                 runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));

           m_spheroids.push_back( STGM::CSpheroid(center,a*s,a,u,s,theta,phi,r,m_spheroids.size()+1,label) );
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
     // set spheroid label
     const char *label = translateChar(asChar(d->label));

     STGM::CSpheroid::direction_type dtype = STGM::CSpheroid::UNIFORM_D;
     if(!strcmp( fname_dir, "rbetaiso" )) {
        dtype=STGM::CSpheroid::BETAISOTROP_D;
     } else if(!strcmp( fname_dir, "rvMisesFisher")) {
        dtype=STGM::CSpheroid::MISES_D;
     }

     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);


     /* shape is always a constant here */
     double shape = asReal(VECTOR_ELT(VECTOR_ELT(d->args,1),0));

     double p1=0,p2=0,
            a=0,b=0,kappa=0,
            theta=0,phi=0,r=0;

     kappa = asReal(VECTOR_ELT(VECTOR_ELT(d->args,2),0));

     p1=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),0));
     if(strcmp(fname, "const" ))
      p2=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),1));

     if(PL>100) {
         Rprintf("Run spheroids  simulation... \n");
         Rprintf("\t size distribution: %s with %f %f \n", fname, p1,p2);
         Rprintf("\t constant shape parameter: %f \n", shape);
         Rprintf("\t directional distribution: %s  with %f \n", fname_dir, kappa);
         Rprintf("\t set label: %s to character: \n",label);
     }

     /* loop over all */
     CVector3d u;
     for (size_t niter=0; niter<num; niter++)
     {
         b = rdist(p1,p2);      /* major semi-axis, for rlnorm possibly switch to perfect simulation */
         a = b * shape;         /* minor semi-axis */
         if(m_stype==CSpheroid::OBLATE)
           std::swap(a,b);

         /* direction */
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

        STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
        m_spheroids.push_back( STGM::CSpheroid(center,a,b,u,shape,theta,phi,r,niter+1,label) );
     }
     PutRNGstate();
}

void STGM::CEllipsoidSystem::IntersectWithPlane(STGM::IntersectorSpheroids &objects, STGM::CPlane &plane)
{
  /// Intersect only objects fully inside the observation window
  //  int i,j;
  //  switch(plane.idx()) {
  //    case 0: i=1; j=2; break; // YZ
  //    case 1: i=0; j=2; break; // XZ
  //    case 2: i=0; j=1; break; // XY
  //  }
  // assume left-down corner is origin of box
  //CWindow win(m_box.m_size[i],m_box.m_size[j]);

  for(size_t i=0; i<m_spheroids.size(); ++i) {
       STGM::IntersectorSpheroid intersector( m_spheroids[i], plane, m_box.m_size);
       if(intersector.FindIntersection()) {
           //if(intersector.getEllipse().isInWindow(win))
             objects.push_back( intersector );
       }
  }

}
