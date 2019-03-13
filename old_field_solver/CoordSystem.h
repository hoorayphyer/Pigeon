#ifndef _COORDSYSTEMB_H_
#define _COORDSYSTEMB_H_

#include <cmath>
#include "Predefs.h"

#define INSTANTIATE_CLASS_WITH_COORDTYPES(ClassName)    \
  template class ClassName<CoordType::CARTESIAN>;       \
  template class ClassName<CoordType::LOG_SPHERICAL>

#define INSTANTIATE_FUNCTION_WITH_COORDTYPES(ret_type, func_name, signature) \
  template ret_type func_name<CoordType::CARTESIAN> signature; \
  template ret_type func_name<CoordType::LOG_SPHERICAL> signature

#define SELECT_COORD(ClassName, type, args...)                          \
  switch (type) {                                                       \
    case CoordType::CARTESIAN:                                          \
      return new ClassName<CoordType::CARTESIAN>;                       \
    case CoordType::LOG_SPHERICAL:                                      \
      return new ClassName<CoordType::LOG_SPHERICAL>;                   \
    default:                                                            \
      return nullptr;                                                   \
  }

// this deals with the case in which FuncName is a template on CoordType
#define SELECT_FUNC_COORD(FuncName, type, args...)     \
  switch (type) {                                      \
  case CoordType::CARTESIAN:                           \
    FuncName<CoordType::CARTESIAN>(args);              \
    break;                                             \
  case CoordType::LOG_SPHERICAL:                       \
    FuncName<CoordType::LOG_SPHERICAL>(args);          \
    break;                                             \
  default:                                             \
    ;                                                  \
  }

using std::sin;
using std::cos;
using std::exp;
using std::acos;


/// Check if the given number is too close to zero. If it is then set
/// it to a small nonzero number
inline double check_eps (double num, double eps = EPS) {
  return ((std::abs(num) < eps) ? std::copysign(eps, num) : num);
}

enum class CoordType : int {
  CARTESIAN = 0,
  LOG_SPHERICAL,
};

struct ScalesCartesian {
  enum { type = static_cast<unsigned int>(CoordType::CARTESIAN) };

  inline double h1(double , double , double = 0.0) const {
    return 1.0;
  }

  inline double h2(double, double, double = 0.0) const {
    return 1.0;
  }

  inline double h3(double, double, double = 0.0) const {
    return 1.0;
  }

  template <int n>
  inline double h(double, double, double = 0.0) const {
    return 1.0;
  }

  // These functions are just decorations and do nothing. They are
  // there to specify a unified interface of this kind of structs
  inline void VectorToCartesian(Vec3<Scalar>& ,
                                   const Vec3<Scalar>& ) const {}

  inline void VectorFromCartesian(Vec3<Scalar>& ,
                                     const Vec3<Scalar>& ) const {}

  inline void PosToCartesian(Vec3<Scalar>& ) const {}

  inline void PosFromCartesian(Vec3<Scalar>& ) const {}
};


struct ScalesLogSpherical {
  enum { type = static_cast<unsigned int>(CoordType::LOG_SPHERICAL) };

  inline double h1(double x, double, double = 0.0) const {
    return exp(x);
  }

  inline double h2(double x, double, double = 0.0) const {
    return exp(x);
  }

  inline double h3(double x, double y, double = 0.0) const {
    return exp(x) * sin(y);
  }

  template <int n>
  inline double h(double r, double theta, double = 0.0) const {
    if (1 == n)
      return exp(r);
    else if (2 == n)
      return exp(r);
    else if (3 == n)
      return exp(r) * check_eps(sin(theta));
    else
      return 1.0;
  }

  // Transform a vector to Cartesian coordinates
  inline void VectorToCartesian(Vec3<Scalar>& v,
                                   const Vec3<Scalar>& pos) const {
    double v1n = v.x, v2n = v.y, v3n = v.z;
    // Now these become vx, vy, vz
    v.x = v1n * sin(pos.y) * cos(pos.z) + v2n * cos(pos.y) * cos(pos.z) -
          v3n * sin(pos.z);
    v.y = v1n * sin(pos.y) * sin(pos.z) + v2n * cos(pos.y) * sin(pos.z) +
          v3n * cos(pos.z);
    v.z = v1n * cos(pos.y) - v2n * sin(pos.y);
  }

  // Transform a vector from Cartesian coordinates
  inline void VectorFromCartesian(Vec3<Scalar>& v,
                                     const Vec3<Scalar>& pos) const {
    double v1n = v.x, v2n = v.y, v3n = v.z;
    // These become vr, vtheta, vphi
    v.x = v1n * sin(pos.y) * cos(pos.z) + v2n * sin(pos.y) * sin(pos.z) +
          v3n * cos(pos.y);
    v.y = v1n * cos(pos.y) * cos(pos.z) + v2n * cos(pos.y) * sin(pos.z) -
          v3n * sin(pos.y);
    v.z = -v1n * sin(pos.z) + v2n * cos(pos.z);
  }

  inline void PosToCartesian(Vec3<Scalar>& pos) const {
    double x1n = pos.x, x2n = pos.y, x3n = pos.z;
    // These become x, y, z
    pos.x = exp(x1n) * sin(x2n) * cos(x3n);
    pos.y = exp(x1n) * sin(x2n) * sin(x3n);
    pos.z = exp(x1n) * cos(x2n);
  }

  inline void PosFromCartesian(Vec3<Scalar>& pos) const {
    double x1n = pos.x, x2n = pos.y, x3n = pos.z;
    // These become r, theta, phi
    pos.x = sqrt(x1n * x1n + x2n * x2n + x3n * x3n);
    pos.y = acos(x3n / pos.x);  // acos(z / r)
    // TODO: Check correctness of this statement!
    pos.z = atan(x2n / x1n) + (x1n < 0) * CONST_PI;  // atan(y / x)
    pos.x = log(pos.x);
  }
};


////////////////////////////////////////////////////////////////////////////////
///  Translating CoordType into the corresponding struct
////////////////////////////////////////////////////////////////////////////////
template <CoordType Coord>
struct CoordToScales {};

template <>
struct CoordToScales<CoordType::CARTESIAN> {
  typedef ScalesCartesian type;
};

template <>
struct CoordToScales<CoordType::LOG_SPHERICAL> {
  typedef ScalesLogSpherical type;
};

#endif  // ----- #ifndef _COORDSYSTEM_H_  -----
