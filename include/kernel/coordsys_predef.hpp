#ifndef  _COORDSYS_PREDEF_HPP_
#define  _COORDSYS_PREDEF_HPP_

namespace knl {
  enum class coordsys : unsigned char
    {
     Cartesian = 0,
     Cylindrical,
     Spherical,
     LogSpherical,
     LogSphericalEV,
    };
}

#endif
