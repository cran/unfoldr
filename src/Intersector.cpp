#include "Intersector.h"

namespace STGM
{

  bool IntersectorSpheroid::TestIntersection ()
  {
    int i=0,j=0,k=0;

    //m_spheroid.CalculateMatrixA(m_iF);
    const CMatrix3d &A = m_spheroid.MatrixA();

    /** Get indices according to the given plane for intersections only by XY,XZ,YZ planes  */
    int l;
    for(l=0; l<3; l++)
        if(m_plane.n[l] == 1 || m_plane.n[l] == -1) break;

    switch (l) {
      case 0: // YZ
        i=1;j=2;k=0;
        break;
      case 1: // XZ
        i=0;j=2;k=1;
        break;
      case 2: // XY
        i=0;j=1;k=2;
        break;
    }

    double d = A[i][i]*A[j][j]-sqr(A[i][j]);
    double s[] = { ( A[i][k]*A[j][j]-A[j][k]*A[i][j] ) /d ,
                  ( A[j][k]*A[i][i]-A[i][k]*A[i][j] ) /d };

    double sum = s[0]*s[0]*A[i][i] + s[0]*s[1]*A[i][j] + s[1]*s[0]*A[j][i] + s[1]*s[1]*A[j][j];

    // translate to xj=0 plane
    CVector3d m(m_spheroid.Center());
    m[k] = m[k] - m_plane.c;
    double tmp = sqr(m[k])*(A[k][k]-sum);

    //condition for intersection
    return ( tmp <= 1);
 }


  bool IntersectorSpheroid::FindIntersection ()
  {
        int i=0,j=0,k=0;
        const CMatrix3d &A = m_spheroid.MatrixA();

        /** Get indices according to the given plane for intersections only by XY,XZ,YZ planes  */
        int l;
        for(l=0; l<3; l++)
            if(m_plane.n[l] == 1 || m_plane.n[l] == -1) break;

        switch (l) {
          case 0: // YZ
            i=1;j=2;k=0;
            break;
          case 1: // XZ
            i=0;j=2;k=1;
            break;
          case 2: // XY
            i=0;j=1;k=2;
            break;
        }

        double d = A[i][i]*A[j][j]-sqr(A[i][j]);
        double s[] = { ( A[i][k]*A[j][j]-A[j][k]*A[i][j] ) /d ,
                      ( A[j][k]*A[i][i]-A[i][k]*A[i][j] ) /d };

        double sum = s[0]*s[0]*A[i][i] + s[0]*s[1]*A[i][j] + s[1]*s[0]*A[j][i] + s[1]*s[1]*A[j][j];

        // translate to xj=0 plane
        CVector3d m(m_spheroid.Center());
        m[k] = m[k] - m_plane.c;
        double tmp = sqr(m[k])*(A[k][k]-sum);

        //condition for intersection
        if( tmp <= 1 )
        {
            CPoint2d m_new;
            CMatrix2d A_new;
            m_new[0]=m[i]+m[k]*s[0];
            m_new[1]=m[j]+m[k]*s[1];
            A_new[0][0]= A[i][i] / (1-tmp);
            A_new[0][1]= A[i][j] / (1-tmp);
            A_new[1][0]= A[j][i] / (1-tmp);
            A_new[1][1]= A[j][j] / (1-tmp);

            /** store ellipse */
            m_ellipse = CEllipse2(A_new,m_new, m_spheroid.Id());

            return true;
       }

      return false;
  }


  bool IntersectorSphere::TestIntersection () {
    double sDist = m_plane.distanceTo(m_sphere.center());
    return fabs(sDist) <= m_sphere.r();
  }

  bool IntersectorSphere::FindIntersection () {
    double sDist = m_plane.distanceTo(m_sphere.center());
    double dist  = fabs(sDist);

    if (dist<=m_sphere.r()) {
      CVector3d center(m_sphere.center());
      center-=sDist*m_plane.n;
      double r = sqrt(fabs(sqr(m_sphere.r())-sqr(dist)));
      m_circle = CCircle3(center,r,m_plane.n,m_sphere.Id());
      return true;
    }
    return false;
  }


} /* namespace STGM */
