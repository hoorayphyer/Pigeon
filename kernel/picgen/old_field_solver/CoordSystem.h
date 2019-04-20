#ifndef _COORDSYSTEMB_H_
#define _COORDSYSTEMB_H_

#include <cmath>
#include "Predefs.h"

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
};

struct ScalesLogSpherical {

  enum { type = static_cast<unsigned int>(CoordType::LOG_SPHERICAL) };

  inline double h1(double x, double, double = 0.0) const {
    return std::exp(x);
  }

  inline double h2(double x, double, double = 0.0) const {
    return std::exp(x);
  }

  inline double h3(double x, double y, double = 0.0) const {
    return std::exp(x) * std::sin(y);
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
