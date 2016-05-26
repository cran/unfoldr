/*
 * IntersectorSpheroidPlane.h
 *
 *  Created on: 31.07.2014
 *      Author: franke
 */

#ifndef INTERSECTOR_H_
#define INTERSECTOR_H_

#include "GeometricPrimitives.h"

namespace STGM
{

  //typedef enum { EMPTY=0 } IntersectionType;
  // Information about the intersection set
  enum IntersectionType  {
    EMPTY=0,          // 0
    NON_EMPTY,        // 1
    POINT,            // 2
    LINE,             // 3
    LINES,            // 4
    CIRCLE,           // 5
    CIRCLE_CAPS,      // 6
    ELLIPSE,          // 7
    ELLIPSE_ARC,      // 8
    ELLIPSE_SEGMENT,  // 9
    ELLIPSE_2D        // 10
  };

  class Intersector
    {
    public:

        virtual ~Intersector () {};

        virtual bool TestIntersection () = 0;
        virtual bool FindIntersection () = 0;

    protected:
        Intersector () {}
    };


  /**
   * @brief IntersectorSpheroidPlane
   *            Only find first intersection with some plane
   */
  class IntersectorSpheroid : public Intersector
  {
  public:
    IntersectorSpheroid(CSpheroid &_spheroid, CPlane &_plane, CPoint3d &_dimensions)
     : m_spheroid(_spheroid), m_plane(_plane), m_size(_dimensions), m_type(EMPTY)
    {
    };

    IntersectorSpheroid(CSpheroid &_spheroid, CPoint3d &_dimensions)
        : m_spheroid(_spheroid), m_size(_dimensions), m_type(EMPTY)
    {
    };

    virtual
    ~IntersectorSpheroid() {}

    bool TestIntersection ();
    bool FindIntersection ();

    void setPlane(CPlane &_plane) { m_plane = _plane; }

    CEllipse2 & getEllipse() { return m_ellipse; }
    const CEllipse2 & getEllipse() const { return m_ellipse; }

    CSpheroid &getSpheroid() { return m_spheroid; }
    const CSpheroid &getSpheroid() const { return m_spheroid; }


    /**
     * @brief Only check if Spheroid intersects a given plane
     *        and store the intersecting plane, translate center coordinate
     *        periodically to the opposite plane
     *
     * @param  plane   [IN]
     * @return boolean [OUT]
     */
    bool operator() (const CPlane &_plane) {
      m_plane = _plane;
      if(TestIntersection()) {
        int l = m_plane.idx();
        m_spheroid.center()[l] +=  m_size[l] * m_plane.n[l];
        return true;
      }
      return false;
    }

    int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
      int interior = 1;
      for(size_t j=0; j<planes.size() ; ++j) {
          if( operator()(planes[j])) {
            interior=0;
            break;
          }
      }
      return interior;
    }

  private:
    CSpheroid m_spheroid;
    CPlane m_plane;

    CPoint3d m_size;
    IntersectionType m_type;

    CEllipse2 m_ellipse;

  };

  typedef std::vector<IntersectorSpheroid> IntersectorSpheroids;

  void digitizeSpheroidIntersections(IntersectorSpheroids &objects, int *w, int nPix, double delta);


  class IntersectorSphere
    : public Intersector
  {
  public:
    IntersectorSphere(CSphere& sphere, CPlane& plane, CPoint3d size) :
       m_sphere(sphere), m_plane(plane), m_size(size), m_type(EMPTY)
    {}

    IntersectorSphere(CSphere &_sphere, CPoint3d &_dimensions) :
       m_sphere(_sphere), m_size(_dimensions), m_type(EMPTY)
    {
    };

    virtual
    ~IntersectorSphere() {}

    bool TestIntersection ();
    bool FindIntersection ();

    void setPlane(CPlane &plane) { m_plane = plane; }

    CCircle3 & getCircle() { return m_circle; }
    const CCircle3 & getCircle() const { return m_circle; }

    CSphere &getSphere() { return m_sphere; }
    const CSphere &getSphere() const { return m_sphere; }

    bool operator() (const CPlane &plane) {
      m_plane = plane;
      if(TestIntersection()) {
        int l = m_plane.idx();
        m_sphere.center()[l] +=  m_size[l] * m_plane.n[l];
        return true;
      }
      return false;
    }

    int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
      int interior = 1;
      for(size_t j=0; j<planes.size() ; ++j) {
          if( operator()(planes[j])) {
            interior=0;
            break;
          }
      }
      return interior;
    }

  private:
    CSphere m_sphere;
    CPlane m_plane;

    CPoint3d m_size;
    IntersectionType m_type;

    CCircle3 m_circle;
  };

  // typedefs
  typedef std::vector<IntersectorSphere> IntersectorSpheres;


  class IntersectorCylinder
      : public Intersector
    {
    public:

      IntersectorCylinder( CCylinder &_cylinder,  CPlane &_plane, CPoint3d &_dimensions)
        : m_cylinder(_cylinder), m_plane(_plane), m_size(_dimensions), m_type(EMPTY), m_side(0)
      {
        setPlaneIdx();
      };


      IntersectorCylinder(  CCylinder &_cylinder, CPoint3d &_dimensions)
        : m_cylinder(_cylinder), m_plane(STGM::CPlane()), m_size(_dimensions), m_type(EMPTY),  m_side(0)
      {
        setPlaneIdx();
      };


      virtual
      ~IntersectorCylinder() {};

      const CCylinder & getCylinder () const { return m_cylinder; };

      CCircle3 & getCircle1 () { return m_circle1; };
      const CCircle3 & getCircle1 () const { return m_circle1; };

      CCircle3 & getCircle2 () { return m_circle2; };
      const CCircle3 & getCircle2 () const { return m_circle2; };

      CEllipse3 & getEllipse () { return m_ellipse; };
      const CEllipse3 & getEllipse () const { return m_ellipse; };

      int getType() const { return m_type; };
      int getSide() const { return m_side; };

      /**
       * @param plane
       * @return
       */
      bool TestIntersection(const CPlane &plane);

      /**
       * @return
       */
      bool TestIntersection () {
        return TestIntersection(m_plane);
      }
      /**
       * @return
       */
      bool FindIntersection ();

      /**
       * @return
       */
      IntersectionType FindIntersectionType();

      inline void setPlaneIdx() {
          switch(m_plane.idx()) {
           case 0: m_i=1; m_j=2; break;
           case 1: m_i=0; m_j=2; break;
           case 2: m_i=0; m_j=1; break;
          }
      }

      /**
       *
       * @param spherecenter
       * @param sDist
       * @return
       */
      CCircle3 GetCircle(STGM::CVector3d &center, double sDist);

      /**
       *
       * @param a
       * @param cosTheta
       * @return
       */
      int FindMajorAxisIntersection(const double a, const double cosTheta);
      /**
       *
       * @param spherecenter
       * @param ipt
       * @return
       */
      double GetEllipseSegment(STGM::CVector3d center, const STGM::CVector3d &ipt);


      /**
       * @brief Only check if Spheroid intersects a given plane
       *        and store the intersecting plane, translate center coordinate
       *        periodically to the opposite plane
       *
       * @param  plane   [IN]
       * @return boolean [OUT]
       */
      bool operator() (const CPlane &plane) {
        if(TestIntersection(plane)) {
          int l = plane.idx();
          m_cylinder.center()[l] +=  m_size[l] * plane.n[l];
          m_cylinder.origin0()[l] = m_cylinder.center()[l] - 0.5*m_cylinder.h()*m_cylinder.u()[l];
          m_cylinder.origin1()[l] = m_cylinder.center()[l] + 0.5*m_cylinder.h()*m_cylinder.u()[l];
          return true;
        }
        return false;
      }

      int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
            int interior = 1;
            for(size_t j=0; j<planes.size() ; ++j) {
                if( operator()(planes[j])) {
                  interior=0;
                  break;
                }
            }
            return interior;
     }

    private:
      CCylinder m_cylinder;
      CPlane m_plane;

      STGM::CPoint3d m_size;
      IntersectionType m_type;

      int m_side, m_i, m_j;
      CCircle3 m_circle1, m_circle2;
      CEllipse3 m_ellipse;

    public:
       STGM::CVector3d ipt0, ipt1;

    };

    /** some type definitions */
    typedef std::vector<IntersectorCylinder> IntersectorCylinders;


    /**
       * @brief Digitize cylinder intersestions
       *
       * @param objects
       * @param w
       * @param nPix
       * @param delta
       */
      void digitizeCylinderIntersections(IntersectorCylinders &objects, int *w, int nPix, double delta);

} /* namespace STGM */

#endif /* INTERSECTOR_H_ */
