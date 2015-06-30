/**
 * @file GeomoetricPrimitives.h
 *
 *  @author: franke
 *  @date  2014-02-11
 */

#ifndef GEOMETRIC_PRIMITIVES_H_
#define GEOMETRIC_PRIMITIVES_H_

#define zero 0

#include <vector>
#include <algorithm>

#include "Vector.h"
#include "Utils.h"

namespace STGM {

  /**
   * Calulate rotation matrix which
   * transforms (0,0,1) into vector v
   */
  CMatrix3d RotationMatrixFrom001(CVector3d v);

  /** some type definitions */
  typedef std::vector<STGM::CPoint2d> PointVector2d;

  /**
   * @brief A base class
   *
   */

  class CWindow;

  class CGeometry
  {
  public:
    CGeometry() {};
    virtual ~CGeometry() {};

  };

   /**
    *  @brief A simple Bounding rectangle aligned to
    *         x and y axis
    */

   class CBoundingRectangle : public CGeometry
   {
     public:
       CBoundingRectangle() :
         m_ymin(0), m_ymax(0), m_xmin(0), m_xmax(0)
       {
       }

       CBoundingRectangle(int ymin, int ymax, int xmin, int  xmax) :
           m_ymin(ymin),m_ymax(ymax), m_xmin(xmin), m_xmax(xmax)
       {
       };

       virtual ~CBoundingRectangle() {};

       int m_ymin, m_ymax;
       int m_xmin, m_xmax;
   };


   /**
    *  @brief Plane
    */
   class CPlane : public CGeometry
   {
   public:

       CPlane () :  n( CVector3d(0,0,1) ), c(0) {}
       virtual ~CPlane () {};

       CPlane ( const CVector3d &_n, const CVector3d &_p) : n(_n)
       {
          c=n.dot(_p);
       }

       CPlane (const CVector3d &_n, const double &_c = 0.0) : n(_n), c(_c) {}

       inline double distanceTo ( const CVector3d &p) const {
         return n.dot(p)-c;
       }

       /**
        * @param p point on plane
        * @return -1,0,1
        */
        inline int whichSide ( const CVector3d &p) const {
          double dist = distanceTo(p);
            if (dist < 0)
                return -1; // p on negative side
            else if (dist > 0)
                return +1; // p on positive side
            else
                return 0;  // p on the plane
        }

        /**
         *
         * @return index of '1'
         */
        inline int idx() const {
           int l=0;
           for(l=0; l<3; l++) if(n[l] == 1 || n[l] == -1) break;
           return l;
        }

        CVector3d n;
        double c;
   };


   /**
    * @brief Circle 3d
    */
   class CCircle3 : public CGeometry
   {
   public:
     CCircle3 () :
       m_center(CVector3d(0,0,0)), m_n(CVector3d(0,0,1)), m_plane(STGM::CPlane()), m_radius(0),m_id(0)
     {
       setPlaneIdx();
     };

     virtual ~CCircle3 (){};

     CCircle3(CVector3d &_center,double &_radius, CVector3d &_n, size_t id)
       : m_center(_center), m_n(_n), m_plane(STGM::CPlane(_n)), m_radius(_radius), m_id(id)
     {
       setPlaneIdx();
     }

     inline void setPlaneIdx() {
        switch(m_plane.idx()) {
          case 0: m_i=1; m_j=2; break; // YZ
          case 1: m_i=0; m_j=2; break; // XZ
          case 2: m_i=0; m_j=1; break; // XY
        }
     }

     void setPlane(CPlane& plane) {
       m_plane =plane;
       m_n = m_plane.n;
       setPlaneIdx();
     }

     size_t Id() { return m_id; }

     CVector3d & center() { return m_center; }
     const CVector3d center() const { return m_center; }

     double & r() { return m_radius; }
     const double r() const { return m_radius; }

     CVector3d& n() { return m_n; }
     const CVector3d n() const { return m_n; }

     inline CPoint2d getMinMax_X() { return CPoint2d(m_center[m_i]-m_radius,m_center[m_i]+m_radius);  }
     inline CPoint2d getMinMax_Y() { return CPoint2d(m_center[m_j]-m_radius,m_center[m_j]+m_radius);  }

     PointVector2d getMinMaxPoints();
     PointVector2d getExtremePointsCircle();

     inline bool isInside(double x, double y)  { return  sqr(x-m_center[m_i])+sqr(y-m_center[m_j])<=sqr(m_radius); };

     CBoundingRectangle & getBoundingRectangle() { return m_br; };

     CVector3d m_center, m_n;
     STGM::CPlane m_plane;
     double m_radius;
     int m_i, m_j;
     CBoundingRectangle m_br;
     size_t m_id;


   };


   class CSphere :  public CGeometry {
      public:

       CSphere(double _x, double _y, double _z, double _r, size_t id=0)
         : type(0), m_id(id), m_center(CVector3d(_x,_y,_z)), m_r(_r)
        {
        }

        CSphere(const CVector3d &_center, const double &_r, const size_t id)
         : type(0), m_id(id), m_center(_center), m_r(_r)
        {
        }

        virtual ~CSphere() {};

        CVector3d& center() { return m_center; }
        const CVector3d& center() const { return m_center; }

        size_t Id() { return m_id; }

        double& r() { return m_r; }
        const double& r() const { return m_r; }

      private:
        size_t type, m_id;
        CVector3d m_center;
        double m_r;

      };

   typedef std::vector<CSphere> Spheres;

  /**
   * @brief Ellipse in 2d
   */

  class CEllipse2 : public CGeometry
  {
  public:

    CEllipse2() :
      m_center(STGM::CPoint2d()),
      m_a(0), m_b(0),
      m_phi(0), id(0), type(10)
    {};

    virtual ~CEllipse2() {};

    /**
     * @brief  Constructor, calculate major and minor axis, angle phi from A
     *
     * @param A_new
     * @param center
     * @param id
     */
    CEllipse2(CMatrix2d &A_new, STGM::CPoint2d &center, size_t id) :
      m_center(center),
      m_A(A_new),
      m_a(0), m_b(0),
      m_phi(0), id(id), type(10)
    {
      int n = 2, err = 0;

      double B[4];
      for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++)
           B[i+n*j] = A_new[i][j];

      double evalf[2] = {0,0};

      /** eigenvalue decomposition */
      real_eval(B,&n,evalf,&err);

      if(!err) {
        double cos_phi = B[1];
        double sin_phi = B[3];

        m_phi = acos(cos_phi);   // angle in the intersecting plane
        if( (cos_phi<0 && sin_phi<0) || (cos_phi<0 && sin_phi>=0)) {
            m_phi = atan(sin_phi/cos_phi)+M_PI;
        } else if(cos_phi>0 && sin_phi<0) {
            m_phi = atan(sin_phi/cos_phi)+2*M_PI;
        }
        m_a = 1/sqrt(evalf[1]);
        m_b = 1/sqrt(evalf[0]);
      } else {
         error(_("CEllipse2(), Error in lapack, Ellipse construction"));
      }


    };

    /**
     * @brief Set dx/dt=0 and dy/dt=0 -> reorder for t values
     *     Helper functions for calculation of extreme points of the ellipse
     *
     * @return Extreme points of the ellipse,
     */
    STGM::CPoint2d getValue() {  return STGM::CPoint2d(atan(-m_b*tan(m_phi)/m_a),atan(m_b/(tan(m_phi)*m_a))); }

    STGM::CPoint2d PointOnEllipse(const double t) {
      return STGM::CPoint2d(m_center[0] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi),
          m_center[1] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi));
    }

    /**
     * @brief Minimum/maximum y coordinates of the ellipse
     *        which is not unique and needs to be sorted
     *
     * @return
     */
    STGM::CPoint2d getMinMax_Y() {
      double t = getValue()[1];
      double y1 = m_center[1] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi);
      double y2 = m_center[1] + m_a*cos(t+M_PI)*sin(m_phi)+m_b*sin(t+M_PI)*cos(m_phi);

      return ( y1<y2 ? STGM::CPoint2d(y1,y2) : STGM::CPoint2d(y2,y1)  );
    }
    /**
     * @brief Minimum/maximum x coordinates of the ellipse
     *       which is not unique and needs to be sorted
     *
     * @return
     */
    STGM::CPoint2d getMinMax_X() {
      double t = getValue()[0];
      double x1 = m_center[0] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi);
      double x2 = m_center[0] + m_a*cos(t+M_PI)*cos(m_phi)-m_b*sin(t+M_PI)*sin(m_phi);

      return ( x1<x2 ? STGM::CPoint2d(x1,x2) : STGM::CPoint2d(x2,x1)  );
    }

   /**
    * @brief Minimum/Maximum coordinates to run through
    *        for digitizing the ellipse intersections
    *
    * @return [min_x,max_x],[min_y,max_y]
    */
    PointVector2d getMinMaxPoints() {
      PointVector2d p;
      p.push_back(getMinMax_X());
      p.push_back(getMinMax_Y());
      return p;
    };


    bool isInside(double x, double y)  {
     double d1 = cos(m_phi)*(x-m_center[0]) + sin(m_phi)*(y-m_center[1]);
     double d2 = sin(m_phi)*(x-m_center[0]) - cos(m_phi)*(y-m_center[1]);

     return (pow(d1 ,2.0)/pow(m_a,2) + pow(d2 ,2.0)/pow(m_b,2)) <= 1;
    }

    /**
     *
     * @return major axis
     */
    double a() const { return m_a; }
    /**
     *
     * @return minor axis
     */
    double b() const { return m_b; }

    /**
     *
     * @return angle major axis to x axis
     */
    double phi() const { return m_phi; }

    /*
    double phi2() {
      double angle = fabs(asin(sin(m_phi)));
      return ( angle<std::numeric_limits<float>::denorm_min() ? 0: angle);
    }
    */

    const STGM::CPoint2d &center()  const { return m_center;}

    /**
     * @return Ellipse matrix
     */
    const CMatrix2d &MatrixA() const { return m_A; }


    inline size_t Id()  const { return id; }

    /**
     * Ellipse2d type
     * @return
     */
    int getType() { return type; }

    /**
    *
    * @return bounding rectangle
    */
    CBoundingRectangle & getBoundingRectangle() { return m_br; }

    /**
    * @brief Does the ellipse intersect the window or
    *        is it fully contained ?
    *
    * @param win Window
    * @return
    */

    /**
    * @brief Extreme points on the ellipse
    *        Not unique, Maximum can return minimum coordinate point
    * @return Point
    */
    STGM::CVector2d getMaxEllipsePoint_X() {
      return STGM::CVector2d(getMinMax_X()[0],getMinMax_Y()[0]);
    }

    STGM::CVector2d getMaxEllipsePoint_Y() {
      return STGM::CVector2d(getMinMax_X()[0],getMinMax_Y()[1]);
    }

    STGM::CVector2d getMinEllipsePoint_X() {
      return STGM::CVector2d(getMinMax_X()[1],getMinMax_Y()[0]);
    }

    STGM::CVector2d getMinEllipsePoint_Y() {
      return STGM::CVector2d(getMinMax_X()[1],getMinMax_Y()[1]);
    }

  private:
    STGM::CPoint2d m_center;
    CMatrix2d m_A;
    double m_a, m_b, m_phi;
    size_t id;
    CBoundingRectangle m_br;
    int type;

  };


   /**
    *   @brief CSpheroid
    */
   class CSpheroid {
     public:

    typedef enum { PROLATE=0, OBLATE=1 } spheroid_type;
     typedef enum {UNIFORM_D=0,BETAISOTROP_D=1, MISES_D=2} direction_type;

     CSpheroid(CVector3d center, double a, double b, CVector3d u,
                 double alpha, double theta, double phi, size_t id )
     :  m_center(center), m_u(u),
        m_a(a), m_b(b),
        m_alpha(alpha),
        m_theta(theta),
        m_phi(phi),
        m_id(id)
      {
        m_R = RotationMatrixFrom001(u);
        m_u.Normalize();
        ComputeMatrixA();
      }

       /**
        * @return rotation matrix
        */
       const CMatrix3d &rotationMatrix() const { return m_R; }

       /**
        *
        * @return
        */
       void ComputeMatrixA();
       const CMatrix3d &MatrixA() const { return m_A; }

       CVector3d &Center() { return m_center; }
       const CVector3d &Center() const { return m_center; }

       double a() const { return m_a; }
       double b() const { return m_b;  }
       double alpha() const { return m_alpha;  }

       const CVector3d& u() const { return m_u; }

       double phi()   const { return m_phi; }
       double theta() const { return m_theta; }

       size_t Id() const { return m_id; }


    private:
      CVector3d m_center, m_u;
      double m_a, m_b;
      double m_alpha,m_theta, m_phi;
      size_t m_id;
      CMatrix3d m_R, m_A, m_invA;

   };

   /**
     *  @brief Box
     */
    class CBox3 : public CGeometry
    {
    public:

      CBox3 ()
       : m_center(STGM::CVector3d(0.5,0.5,0.5)),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(STGM::CPoint3d(1,1,1))
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5;
        m_extent[1] = 0.5;
        m_extent[2] = 0.5;

        ConstructBoundingPlanes();
      }

      CBox3 (double xrange[2], double yrange[2], double zrange[2])
       :  m_u(STGM::CVector3d(1.0,0.0,0.0)),
          m_v(STGM::CVector3d(0.0,1.0,0.0)),
          m_w(STGM::CVector3d(0.0,0.0,1.0)),
          m_size(STGM::CPoint3d(fabs(xrange[1]-xrange[0]),fabs(yrange[1]-yrange[0]),fabs(zrange[1]-zrange[0])))
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5*m_size[0];
        m_extent[1] = 0.5*m_size[1];
        m_extent[2] = 0.5*m_size[2];

        m_center[0] = xrange[1]-m_extent[0];
        m_center[1] = yrange[1]-m_extent[1];
        m_center[2] = zrange[1]-m_extent[2];

        ConstructBoundingPlanes();
      }

      CBox3 (STGM::CVector3d center, STGM::CPoint3d size)
       : m_center(center),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(size)
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5*m_size[0];
        m_extent[1] = 0.5*m_size[1];
        m_extent[2] = 0.5*m_size[2];

        ConstructBoundingPlanes();
      }

      /** @brief Constructor Box from left lower point (0,0,0)
       *
       * @param a  extent in x in first direction
       * @param b  extent in y in second direction
       * @param c  extent in y in third direction
       */

      CBox3 (double a, double b, double c)
      :  m_center(STGM::CVector3d(0.5*a,0.5*b,0.5*c)),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(STGM::CPoint3d(a,b,c)) /// length in each direction
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        /// half lengths
        m_extent[0] = 0.5*a;
        m_extent[1] = 0.5*b;
        m_extent[2] = 0.5*c;

        ConstructBoundingPlanes();
      }

      virtual ~CBox3 () {};

      // public members
      void setExtent(double a, double b, double c);

      double volume() const { return  m_size[0]*m_size[1]*m_size[2]; };

      void ConstructBoundingPlanes();

      double PointInBox3(STGM::CVector3d &point);

      const std::vector<CPlane> &getPlanes() const { return m_planes; };

      /**
       * @brief Get planes perpendicular to the argument plane
       *
       * @param plane
       * @return planes
       */
      const std::vector<CPlane> getPlanes( CPlane &plane) const {
        std::vector<CPlane> planes;
        int l = plane.idx();
        for(int k=0; k<4; ++k) {
            if(k!=2*l && k!=2*l+1)
              planes.push_back(m_planes[k]);
        }
           return planes;
      }

      STGM::CVector3d m_center;
      STGM::CVector3d m_u, m_v, m_w, *m_axis[3];
      double m_extent[3];

      STGM::CPoint3d m_size;
      std::vector<CPlane> m_planes;
    };


   typedef std::vector<CSpheroid> Spheroids;
   typedef std::vector<CEllipse2> Ellipses;


} /* namespace STGM */

#endif /* GEOMETRIC_PRIMITIVES_H_ */
