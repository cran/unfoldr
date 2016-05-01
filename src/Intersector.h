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

typedef enum { EMPTY=0 } IntersectionType;

  /**
   * @brief IntersectorSpheroidPlane
   *            Only find first intersection with some plane
   */
  class IntersectorSpheroid
  {
  public:
    IntersectorSpheroid(CSpheroid &_spheroid, CPlane &_plane, CPoint3d &_dimensions)
     : m_spheroid(_spheroid), m_plane(_plane), m_size(_dimensions)
    {
      m_type = EMPTY;
    };

    IntersectorSpheroid(CSpheroid &_spheroid, CPoint3d &_dimensions)
        : m_spheroid(_spheroid), m_size(_dimensions)
    {
      m_type = EMPTY;
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
        m_spheroid.Center()[l] +=  m_size[l] * m_plane.n[l];
        return true;
      }
      return false;
    }

  private:
    IntersectionType m_type;
    CSpheroid m_spheroid;
    CPlane m_plane;
    CPoint3d m_size;
    CEllipse2 m_ellipse;

  };

  typedef std::vector<IntersectorSpheroid> IntersectorSpheroidPlaneVec;

  void digitizeSpheroidIntersections(IntersectorSpheroidPlaneVec &objects, int *w, int nPix, double delta);


  class IntersectorSphere
  {
  public:
    IntersectorSphere(CSphere& sphere, CPlane& plane, CPoint3d size) :
       m_sphere(sphere), m_plane(plane), m_size(size)
    {}

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

  private:

    CSphere m_sphere;
    CPlane m_plane;
    CPoint3d m_size;
    CCircle3 m_circle;
  };

  typedef std::vector<IntersectorSphere> IntersectorSpherePlaneVec;

} /* namespace STGM */

#endif /* INTERSECTOR_H_ */
