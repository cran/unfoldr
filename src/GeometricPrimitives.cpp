/** @file  GeometricPrimitives.cpp
 *
 *  @date 02-14-2014
 *  @uthor: franke
 */
#define ZERO_TOL 1E-6

#define NDIM 3
#define DISTANCE_VEC(P1,P2,R)                      \
do {                                               \
  int _i;                                          \
  for (_i = 0; _i < NDIM; _i++)                    \
    (R)[_i] = (P1)[_i]-(P2)[_i];                   \
} while(0)


#include <R_ext/Lapack.h>
#include "GeometricPrimitives.h"

static inline double sign(double a,double b) { return a = fabs(a),(b<0)?-a:a; }

namespace STGM {
  typedef const double *array_t;
  typedef const double (*matrix_t)[3];

  void real_eval(double *a, int *n, double *evalf, int *err) {
              // size(evalf): n
              // result: evec = a

              int lda = *n,  lwork = 3*lda-1;
              double *work = Calloc(lwork, double);
              F77_NAME(dsyev)("V","U", &lda, a, &lda, evalf, work, &lwork, err);
              Free(work);

  }

  double CWindow::PointInWindow(STGM::CVector2d point)  {
     STGM::CVector2d diff(point);
     diff -= m_center;

     double sqrDistance = 0, delta=0;
     STGM::CVector2d closest;

     for (int i = 0; i < 2; ++i)
     {
         closest[i] = diff.dot( *(m_axis[i]));
         if (closest[i] < -m_extent[i])
         {
             delta = closest[i] + m_extent[i];
             sqrDistance += delta*delta;
             closest[i] = -m_extent[i];
         }
         else if (closest[i] > m_extent[i])
         {
             delta = closest[i] - m_extent[i];
             sqrDistance += delta*delta;
             closest[i] = m_extent[i];
         }
     }

     STGM::CVector2d mClosestPoint = m_center;
     for (int i = 0; i < 2; ++i)
       mClosestPoint += closest[i]* *(m_axis[i]);

     return sqrDistance;
  }


  CVector3d PerpendicularVector(const CVector3d &a) {
      if (fabs(a[2]) > fabs(a[0]))
        return CVector3d(0.0, -a[2], a[1]);
      else
        return CVector3d(-a[1], a[0], 0.0);
  }

  CMatrix3d RotationMatrixFrom001(CVector3d v) {
    v.Normalize();

    CVector3d o1(PerpendicularVector(v));
    o1.Normalize();

    CVector3d o2(cross(v,o1));
    o2.Normalize();

    CMatrix3d r;

    r[0][0] = o1[0];
    r[0][1] = o1[1];
    r[0][2] = o1[2];

    r[1][0] = o2[0];
    r[1][1] = o2[1];
    r[1][2] = o2[2];

    r[2][0] = v[0];
    r[2][1] = v[1];
    r[2][2] = v[2];

    return r;
  }

  double shortestDistance(const double *r12,  const double *u1, const double *u2, const double &lh1,
                            const double &lh2, double &xla, double &xmu)
  {
    double  rr =  r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2],
            ru1 = r12[0]*u1[0]+r12[1]*u1[1]+r12[2]*u1[2],
            ru2 = r12[0]*u2[0]+r12[1]*u2[1]+r12[2]*u2[2],
            u1u2 = u1[0]*u2[0]+u1[1]*u2[1]+u1[2]*u2[2],
            cc = 1.0-SQR(u1u2);

    // Checking whether the rods are or not parallel:
    // The original code is modified to have symmetry:

     if(cc<ZERO_TOL) {
      if(ru1!=0 && ru2!=0) {
          xla= ru1/2;
          xmu= -ru2/2;
      } else {
          xla=0;
          xmu=0;
          return rr;
      }
     } else {
      // Step 1
      xla= (ru1-u1u2*ru2)/cc;
      xmu= (-ru2+u1u2*ru1)/cc;
     }

    // Step 2
    if( fabs(xla)>lh1 || fabs(xmu)>lh2 ) {
    // Step 3 - 7
      if(fabs(xla)-lh1>fabs(xmu)-lh2) {
       xla= sign(lh1,xla);
       xmu= xla*u1u2-ru2;
       if( fabs(xmu)>lh2 ) xmu= sign(lh2,xmu);
      }
      else {
       xmu= sign(lh2,xmu);
       xla= xmu*u1u2+ru1;
       if( fabs(xla)>lh1 ) xla= sign(lh1,xla);
      }
     }
     // Step 8
      return rr+SQR(xla)+SQR(xmu) - 2*xla*xmu*u1u2 + 2*xmu*ru2 - 2*xla*ru1;
  }


  double CSpheroid::spheroidDistanceAsCylinder(CSpheroid &sp) const {
      double d=0, xla=0,xmu=0;
      double lh1=1.0e-7,lh2=1.0e-7;

      STGM::CPoint3d w(0,0,0);

      //sp.setCrackType(crack);
      DISTANCE_VEC(m_center.ptr(),sp.Center().ptr(),w.ptr());

      if(m_crack && sp.m_crack)  {  // E_i, E_j
         lh1 = m_b-m_a;
         lh2 = sp.b()-sp.a();
         d=shortestDistance(w.ptr(),m_u.ptr(),sp.u().ptr(),lh1,lh2,xla,xmu);
         return MAX(sqrt(d) - m_a - sp.a(),0.0);

      } else if(!m_crack && sp.m_crack)  {  // D_i, E_j
         lh2 = sp.b()-sp.a();
         d=shortestDistance(w.ptr(),m_u.ptr(),sp.u().ptr(),lh1,lh2,xla,xmu);
         return MAX(sqrt(d) - sp.a(),0.0);

      } else if(m_crack && !sp.m_crack)  {  // E_i, D_j
          lh1 = m_b-m_a;
          d=shortestDistance(w.ptr(),m_u.ptr(),sp.u().ptr(),lh1,lh2,xla,xmu);
          return MAX(sqrt(d) - m_a,0.0);

      } else  {//if(!m_crack && !sp.m_crack)  // D_i, D_j
          d=sqrt(shortestDistance(w.ptr(),m_u.ptr(),sp.u().ptr(),lh1,lh2,xla,xmu));
          return MAX(d,0.0);
      }

    }


  void CSpheroid::ComputeMatrixA() {
     m_A[0][0] = m_A[1][1] = 1.0 / SQR(m_a);
     m_A[2][2] = 1.0 / SQR(m_b);

     CMatrix3d R = RotationMatrixFrom001(m_u);

     m_A = m_A * R;
     R.Transpose();
     m_A = R * m_A;
   }

   PointVector2d CCircle3::getMinMaxPoints()  {
     PointVector2d p;
     p.push_back(getMinMax_X());
     p.push_back(getMinMax_Y());
     return p;
   }

   PointVector2d CCircle3::getExtremePointsCircle() {
     PointVector2d p;
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[0]));
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[1]));
     p.push_back(CPoint2d(getMinMax_X()[0],m_center[m_j]));
     p.push_back(CPoint2d(getMinMax_X()[1],m_center[m_j]));
     return p;
   }

   /*
   bool CCircle3::isInWindow(STGM::CWindow &win) {
    std::vector<CPoint2d> p = getMinMaxPoints();
    if( (win.PointInWindow( CVector2d(p[0][0],p[1][0]) ) == 0) &&
       (win.PointInWindow( CVector2d(p[0][0],p[1][1]) ) == 0) &&
       (win.PointInWindow( CVector2d(p[0][1],p[1][0]) ) == 0) &&
       (win.PointInWindow( CVector2d(p[0][1],p[1][1]) ) == 0))
    {
     return true;
    }
   return false;
  }
  */


  void CBox3::setExtent(double a, double b, double c) {
     m_extent[0] = 0.5*a;
     m_extent[1] = 0.5*b;
     m_extent[2] = 0.5*c;

     m_size[0] = a;
     m_size[1] = b;
     m_size[2] = c;
  }

  void CBox3::ConstructBoundingPlanes()
  {
    /** sorted order of normal vectors (planes) -  do not change! */
    m_planes.push_back(CPlane(CVector3d(1,0,0),0)); // left
    m_planes.push_back(CPlane(CVector3d(-1,0,0),m_size[0])); // right
    m_planes.push_back(CPlane(CVector3d(0,1,0),0)); // front
    m_planes.push_back(CPlane(CVector3d(0,-1,0),m_size[1])); // behind
    m_planes.push_back(CPlane(CVector3d(0,0,1),0)); // bottom
    m_planes.push_back(CPlane(CVector3d(0,0,-1),m_size[2])); // top
  }

  void CBox3::ConstructBoxLateralPlanes()
  {
      /** sorted order of normal vectors (planes) -  do not change! */
      m_lateral_planes.push_back(CPlane(STGM::CVector3d(1,0,0),0)); // left
      m_lateral_planes.push_back(CPlane(STGM::CVector3d(-1,0,0),m_size[0])); // right
      m_lateral_planes.push_back(CPlane(STGM::CVector3d(0,1,0),0)); // front
      m_lateral_planes.push_back(CPlane(STGM::CVector3d(0,-1,0),m_size[1])); // behind
  }

  STGM::CEllipse2 CSpheroid::spheroidCrackProjection() const {
    if(m_crack)
      return spheroidProjection();
    else
      return spheroidCircleProjection();
  }


  STGM::CEllipse2 CSpheroid::spheroidProjection() const
  {
    array_t  m = m_center.ptr();

    STGM::CMatrix2d A_new;
    A_new[0][0]= m_A[0][0];
    A_new[0][1]= m_A[0][1];
    A_new[1][0]= m_A[1][0];
    A_new[1][1]= m_A[1][1];

    STGM::CPoint2d center(m[0],m[1]);
    return STGM::CEllipse2(A_new, center, m_id, 0.5*M_PI);
  }

  STGM::CEllipse2 CSpheroid::spheroidCircleProjection() const
  {
        STGM::CPoint2d ctr(m_center[0],m_center[1]);
        STGM::CVector3d n(0,0,1), u(m_u);

        STGM::CVector3d z(cos(m_phi)*u[2],
                          sin(m_phi)*u[2],
                          sin(-m_phi)*u[1]-cos(-m_phi)*u[0]);

        STGM::CVector3d pz(z),               // minor
                        pv(cross(m_u,z));    // major

        pz.Normalize();
        pv.Normalize();

        pz *= m_a;
        pz += m_center;
        pv *= m_a;
        pv += m_center;

        STGM::CPoint2d axis1(pz[0]-ctr[0],pz[1]-ctr[1]);
        STGM::CPoint2d axis2(pv[0]-ctr[0],pv[1]-ctr[1]);

        double zlen = axis1.Length();
        double vlen = axis2.Length();

        axis1.Normalize();
        axis2.Normalize();

        return STGM::CEllipse2(ctr,axis2,axis1,zlen,vlen,m_id);
  }

} /* STGM */
