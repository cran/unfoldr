#ifndef SRC_SIMSPHERE_H_
#define SRC_SIMSPHERE_H_

#include "Intersector.h"

#ifdef __cplusplus
extern "C" {
#endif

  SEXP SphereSystem(SEXP R_param, SEXP R_cond);
  SEXP GetSphereSystem(SEXP ext);
  SEXP SetupSphereSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond);
  SEXP IntersectSphereSystem(SEXP ext, SEXP R_n, SEXP R_z);
  SEXP FinalizeSphereSystem(SEXP ext);

#ifdef __cplusplus
}
#endif


namespace STGM {

class CBoolSphereSystem {
public:

  CBoolSphereSystem(CBox3 &box, double lam) :
    m_box(box),m_lam(lam), num(0)
  {}

    void simSphereSys(R_Calldata d);
    size_t getNumSpheres() const { return num; }
    STGM::Spheres &refObjects()  { return m_spheres; }
    const STGM::Spheres &refObjects() const { return m_spheres; }

    void IntersectWithPlane(IntersectorSpherePlaneVec &objects, STGM::CPlane &plane);

    template<typename F> void simSpheres(F f);



private:
  CBox3 m_box;
  double m_lam;
  size_t num;
  Spheres m_spheres;

};

}


#endif /* SRC_SIMSPHERE_H_ */
