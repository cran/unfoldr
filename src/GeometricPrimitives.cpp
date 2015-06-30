/** @file  GeometricPrimitives.cpp
 *
 *  @date 02-14-2014
 *  @uthor: franke
 */

#include "GeometricPrimitives.h"

namespace STGM {

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


  void CSpheroid::ComputeMatrixA() {
     m_A[0][0] = m_A[1][1] = 1.0 / sqr(m_a);
     m_A[2][2] = 1.0 / sqr(m_b);

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
     //Rprintf("getExtremePointsCircle...\n");
     PointVector2d p;
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[0]));
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[1]));
     p.push_back(CPoint2d(getMinMax_X()[0],m_center[m_j]));
     p.push_back(CPoint2d(getMinMax_X()[1],m_center[m_j]));
     return p;
   }


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
    m_planes.push_back(CPlane(CVector3d(1,0,0),zero)); // left
    m_planes.push_back(CPlane(CVector3d(-1,0,0),m_size[0])); // right
    m_planes.push_back(CPlane(CVector3d(0,1,0),zero)); // front
    m_planes.push_back(CPlane(CVector3d(0,-1,0),m_size[1])); // behind
    m_planes.push_back(CPlane(CVector3d(0,0,1),zero)); // bottom
    m_planes.push_back(CPlane(CVector3d(0,0,-1),m_size[2])); // top
  }


} /* STGM */
