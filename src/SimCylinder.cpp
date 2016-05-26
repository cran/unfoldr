/**
 * @file SimCylinder.h
 *
 *  @date: 09.05.2016
 *  @author: M. Baaske
 */

#include <vector>

#include "directions.h"
#include "SimCylinder.h"

static SEXP cylinder_type_tag;
static int PL = 0;

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_i+DIM*_j] = (M)[_i][_j];         \
} while(0)

#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_i+DIM*_j];         \
} while(0)

SEXP convert_R_Cylinder( STGM::CCylinder &cyl, STGM::LateralPlanes &planes , STGM::CBox3 &box);
SEXP convert_R_Cylinders( STGM::Cylinders &cyl, STGM::CBox3 &box);
STGM::CCylinder convert_C_Cylinder(SEXP R_cyl);

void _cylp_finalizer(SEXP ext)
{
    checkPtr(ext, cylinder_type_tag);
    STGM::CCylinderSystem *sptr = static_cast<STGM::CCylinderSystem *>(R_ExternalPtrAddr(ext));
    delete sptr;
    R_ClearExternalPtr(ext);
}


SEXP FinalizeCylinderSystem(SEXP ext) {
  _cylp_finalizer(ext);
  return R_NilValue;
}

SEXP CreateExternalCylinderPointer(STGM::CCylinderSystem *sp) {
  cylinder_type_tag = install("CylinderSystem_TAG");
  SEXP Rextp=R_NilValue;
  PROTECT(Rextp=R_MakeExternalPtr((void*) sp, cylinder_type_tag, R_NilValue));
  R_RegisterCFinalizerEx(Rextp, (R_CFinalizer_t) _cylp_finalizer, TRUE);

  UNPROTECT(1);
  return Rextp;
}

SEXP UpdateIntersectionsCylinder(SEXP R_cylinder, SEXP R_cluster_ids, SEXP R_box) {
    int nProtected=0;
    SEXP R_ret, R_clust_id;

    SEXP R_BoxX, R_BoxY, R_BoxZ;
    PROTECT( R_BoxX = AS_NUMERIC( getListElement( R_box, "xrange" ) ) ); ++nProtected;
    PROTECT( R_BoxY = AS_NUMERIC( getListElement( R_box, "yrange" ) ) ); ++nProtected;
    PROTECT( R_BoxZ = AS_NUMERIC( getListElement( R_box, "zrange" ) ) ); ++nProtected;

    STGM::CBox3 box(REAL(R_BoxX)[1],REAL(R_BoxY)[1],REAL(R_BoxZ)[1]);
    const std::vector<STGM::CPlane> &planes = box.getPlanes();

    PROTECT(R_clust_id = getListElement(R_cluster_ids,"id")); ++nProtected;
    PROTECT(R_ret = allocVector(INTSXP,length(R_clust_id))); ++nProtected;

    for(int k=0;k<length(R_clust_id);k++) {
        STGM::CCylinder sp = convert_C_Cylinder(VECTOR_ELT(R_cylinder,k));
        STGM::IntersectorCylinder intersector(sp , box.m_size );
        for(size_t j=0; j<planes.size() ; ++j) {
            if( intersector(planes[j])) {
              sp.interior()=0;
              break;
            }
        }
        INTEGER(R_ret)[k]=sp.interior();
    }
    UNPROTECT(nProtected);
    return R_ret;
  }



STGM::CCylinderSystem * InitCylinderSystem(SEXP R_param, SEXP R_cond) {
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  // set print level
  PL = asInteger(getListElement( R_cond,"pl"));
  double lam = asReal(AS_NUMERIC(getListElement( R_param, "lam")));

  /* set up cylinder system */
  STGM::CBox3 box(boxX,boxY,boxZ);
  STGM::CVector3d maxis( NUMERIC_POINTER( getListElement( R_cond, "mu")) );

  STGM::CCylinderSystem *sp = (STGM::CCylinderSystem *)Calloc(1,STGM::CCylinderSystem);

  try {
      new(sp)STGM::CCylinderSystem(box,lam,maxis);
  } catch(...) {
      error(_("InitCylinderSystem(): Allocation error."));
  }

  UNPROTECT(1);
  return sp;
}

SEXP CylinderSystem(SEXP R_param, SEXP R_cond)
{
  STGM::CCylinderSystem *sp = InitCylinderSystem(R_param,R_cond);
  R_Calldata call_data = getRCallParam(R_param,R_cond);

  if(TYPEOF(call_data->fname) != VECSXP) {
        sp->simSysJoint(call_data);
    } else {
        if(!strcmp(GET_NAME(call_data,0), "rbinorm") ) {
            sp->simBivariate(call_data);
        } else {
            sp->simCylinderSys(call_data);
        }
  }
  deleteRCall(call_data);

  if(PL>100) {
    Rprintf("Simulated %d spheroids: %p \n",sp->refObjects().size(), sp);
  }

  SEXP R_cylinders = R_NilValue;
  if(PL>100) {
      STGM::Cylinders &cylinders = sp->refObjects();
      Rprintf("Convert... \n");
      PROTECT(R_cylinders = convert_R_Cylinders( cylinders, sp->box() ));
  } else {
      PROTECT(R_cylinders = allocVector(VECSXP,0));  /* return empty list */
  }

  SEXP Rextp=R_NilValue;
  PROTECT(Rextp = CreateExternalCylinderPointer(sp));

  setAttrib(R_cylinders, install("eptr"), Rextp);
  classgets(R_cylinders, mkString("cylinder"));

  UNPROTECT(2);
  return R_cylinders;
}

void STGM::CCylinderSystem::simSysJoint(R_Calldata d) {
     GetRNGstate();

     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
       num = rpois(m_box.volume()*m_lam);
       ++nTry;
     }
     m_cylinders.reserve(num);
     // set cylinder label
     const char *label = translateChar(asChar(d->label));
     // set box constants
     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double h=0,radius=0,theta=0, phi=0,r=0;
     double *v;

     SEXP Reval=R_NilValue;
     CVector3d u;

     int err = 0;
     for (size_t niter=0; niter<num; niter++) {
         Reval = R_tryEval(d->call,d->rho,&err);
         if(!err) {
            h=asReal(getListElement(Reval,"h"));
            r=asReal(getListElement(Reval,"radius"));
            theta=asReal(getListElement(Reval,"theta"));
            phi=asReal(getListElement(Reval,"phi"));

            v=REAL(getListElement(Reval,"u"));
            u[0]=v[0]; u[1]=v[1];  u[2]=v[2];

            STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
            m_cylinders.push_back(STGM::CCylinder(center,u,h,r,theta,phi,radius,niter+1,label) );
         } else
           error(_("simSysJoint(): R try error in user defined distribution function."));
     }
     PutRNGstate();
}

void STGM::CCylinderSystem::simBivariate(R_Calldata d) {
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
  m_cylinders.reserve(num);

  CVector3d u;
  double x=0,y=0,h=0,s=1,phi=0,theta=0;

  int k=0;
  double r=0;
  for (size_t niter=0; niter<num; niter++)
  {
      /* sample major semi-axis a, shorter semi-axis is: c=a*s */
      rbinorm(mx,sdx,my,sdy,rho,x,y);
      s=1.0/(1.0+exp(-y));
      h=exp(x);

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

       m_cylinders.push_back(STGM::CCylinder(center,u,h,h*s,theta,phi,r,niter+1,label) );
   }
   PutRNGstate();

}

void STGM::CCylinderSystem::simCylinderSys(R_Calldata d) {
  GetRNGstate();

  int nTry=0;
  while(num==0 && nTry<MAX_ITER) {
     num = rpois(m_box.volume()*m_lam);
     ++nTry;
  }
  m_cylinders.reserve(num);

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

  Rboolean isRadius = FALSE;
  const char *fname_shape = GET_NAME(d,1);
  if(!strcmp(fname_shape, "radius" )) {
    isRadius = TRUE;
  }

  const char *fname_dir = GET_NAME(d,2);
  const char *label = translateChar(asChar(d->label));

  STGM::CCylinder::direction_type dtype = STGM::CCylinder::UNIFORM_D;
  if(!strcmp( fname_dir, "rbetaiso" )) {
     dtype=STGM::CCylinder::BETAISOTROP_D;
  } else if(!strcmp( fname_dir, "rvMisesFisher")) {
     dtype=STGM::CCylinder::MISES_D;
  }

  double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
         m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
         m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);


  /* shape is always a constant here */
  double shape = asReal(VECTOR_ELT(VECTOR_ELT(d->args,1),0));

  double p1=0, p2=0,
         h=0, radius=0,
         kappa=0, theta=0, phi=0, r=0;

  // orientation iso parameter
  kappa = asReal(VECTOR_ELT(VECTOR_ELT(d->args,2),0));
  // distribution parameters
  p1=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),0));
  if(strcmp(fname, "const" ))
   p2=asReal(VECTOR_ELT(VECTOR_ELT( d->args, 0),1));

  if(PL>100) {
      Rprintf("Run cylinder simulation... \n");
      Rprintf("\t size distribution: %s with %f %f \n", fname, p1,p2);

      if(isRadius)
        Rprintf("\t constant radius: %f \n", shape);
      else
        Rprintf("\t constant shape parameter: %f \n", shape);

      Rprintf("\t directional distribution: %s  with %f \n", fname_dir, kappa);
      Rprintf("\t set label: %s to character: \n",label);
  }

  /* loop over all */
  CVector3d u;
  for (size_t niter=0; niter<num; niter++)
  {
      h = rdist(p1,p2);                                /* cylinder axis, for rlnorm possibly switch to perfect simulation */
      r = (isRadius ? shape : h * shape);         /* either fixed radius or radius as portion of height */

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
     m_cylinders.push_back( STGM::CCylinder(center,u,h,r,theta,phi,radius,niter+1,label) );
  }
  PutRNGstate();
}

STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  double radius = 0;
  const char *label = "N";

  if(!isNull(getAttrib(R_cyl, install("label"))))
    label = translateChar(asChar(getAttrib(R_cyl, install("label"))));
  if(!isNull(getAttrib(R_cyl, install("interior"))))
    interior = asLogical(getAttrib(R_cyl, install("interior")));
  if(!isNull(getAttrib(R_cyl, install("radius"))))
    radius = asReal(getAttrib(R_cyl, install("radius")));

  if(!strcmp(label,"F")) {
      SEXP R_ctr;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      STGM::CVector3d ctr(REAL(R_ctr)),u(0,0,1);

      UNPROTECT(1);
      return STGM::CCylinder(ctr,u,0,asReal(getListElement(R_cyl, "r")),0,0,
                 radius, asInteger(getListElement(R_cyl, "id")), label, interior);

  } else {
      SEXP R_ctr, R_u, R_angles;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_cyl, "u")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_cyl, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)),u(REAL(R_u));
      UNPROTECT(3);

      return STGM::CCylinder(ctr,u, asReal(getListElement(R_cyl, "length")),
                asReal(getListElement(R_cyl, "r")), REAL(R_angles)[0],REAL(R_angles)[1],
                  radius, asInteger(getListElement(R_cyl, "id")), label, interior);
  }

}


SEXP convert_R_Cylinder( STGM::CCylinder &cyl, STGM::LateralPlanes &planes, STGM::CBox3 &box) {
  int ncomps=9, dim=3;
  SEXP names, R_Cyl = R_NilValue;
  SEXP R_u, R_center,  R_origin0, R_origin1, R_angles, R_rotM;

  PROTECT(R_Cyl = allocVector(VECSXP,ncomps));
  PROTECT(R_u = allocVector(REALSXP, 3));
  PROTECT(R_center = allocVector(REALSXP, 3));
  PROTECT(R_origin0 = allocVector(REALSXP, 3));
  PROTECT(R_origin1 = allocVector(REALSXP, 3));
  PROTECT(R_angles = allocVector(REALSXP,2));
  PROTECT(R_rotM = allocMatrix(REALSXP,3,3));

  REAL(R_u)[0]=cyl.u()[0];
  REAL(R_u)[1]=cyl.u()[1];
  REAL(R_u)[2]=cyl.u()[2];

  REAL(R_center)[0]=cyl.center()[0];
  REAL(R_center)[1]=cyl.center()[1];
  REAL(R_center)[2]=cyl.center()[2];

  REAL(R_origin0)[0] = cyl.origin0()[0];
  REAL(R_origin0)[1] = cyl.origin0()[1];
  REAL(R_origin0)[2] = cyl.origin0()[2];

  REAL(R_origin1)[0] = cyl.origin1()[0];
  REAL(R_origin1)[1] = cyl.origin1()[1];
  REAL(R_origin1)[2] = cyl.origin1()[2];

  REAL(R_angles)[0]=cyl.theta();
  REAL(R_angles)[1]=cyl.phi();

  //for (size_t i = 0; i < 3; i++)
  //  for (size_t j = 0; j < 3; j++)
  //    REAL(R_rotM)[i + 3 *j] = (cyl.rotationMatrix())[i][j];

  const STGM::CMatrix3d &M = cyl.rotationMatrix();
  COPY_C2R_MATRIX(M,R_rotM,dim);

  // check intersection
  STGM::IntersectorCylinder intersector(cyl, box.m_size );
  Rboolean interior = (Rboolean) TRUE;
  for(size_t j=0; j<planes.size() ; ++j) {
       if( intersector(planes[j])) {
         interior = (Rboolean) FALSE;
         break;
       }
  }
  STGM::PointVector2d p;
  double area = cyl.delamProjection(p,20);

  SET_VECTOR_ELT(R_Cyl,0,ScalarInteger( (int) cyl.Id() ));
  SET_VECTOR_ELT(R_Cyl,1,R_center);
  SET_VECTOR_ELT(R_Cyl,2,R_origin0);
  SET_VECTOR_ELT(R_Cyl,3,R_origin1);
  SET_VECTOR_ELT(R_Cyl,4,ScalarReal(cyl.h())) ;
  SET_VECTOR_ELT(R_Cyl,5,R_u);
  SET_VECTOR_ELT(R_Cyl,6,ScalarReal(cyl.r()));
  SET_VECTOR_ELT(R_Cyl,7,R_angles);
  SET_VECTOR_ELT(R_Cyl,8,R_rotM);

  PROTECT(names = allocVector(STRSXP, ncomps));
  SET_STRING_ELT(names, 0, mkChar("id"));
  SET_STRING_ELT(names, 1, mkChar("center"));
  SET_STRING_ELT(names, 2, mkChar("origin0"));
  SET_STRING_ELT(names, 3, mkChar("origin1"));
  SET_STRING_ELT(names, 4, mkChar("length"));
  SET_STRING_ELT(names, 5, mkChar("u"));
  SET_STRING_ELT(names, 6, mkChar("r"));
  SET_STRING_ELT(names, 7, mkChar("angles"));
  SET_STRING_ELT(names, 8, mkChar("rotM"));

  setAttrib(R_Cyl, R_NamesSymbol, names);
  setAttrib(R_Cyl, install("radius"), ScalarReal(cyl.radius()));
  setAttrib(R_Cyl, install("label"), mkString(cyl.label()) );
  setAttrib(R_Cyl, install("interior"), ScalarLogical(interior));
  // here always 'delam' projection area
  setAttrib(R_Cyl, install("area"), ScalarReal(area));

  UNPROTECT(8);
  return R_Cyl;
}

SEXP convert_R_Cylinders( STGM::Cylinders &cyls, STGM::CBox3 &box) {
  Rprintf("Convert...%d  cylinders.\n", cyls.size());

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, cyls.size()) );
  // get lateral bounding planes
  STGM::LateralPlanes &planes = box.getLateralPlanes();

  for(size_t k=0;k<cyls.size();++k)
    SET_VECTOR_ELT(R_resultlist,k,convert_R_Cylinder(cyls[k],planes,box));

  UNPROTECT(1);
  return(R_resultlist);
}

