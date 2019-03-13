#ifndef _TYPES_H_
#define _TYPES_H_

#include "Predefs.h"
#include <cstdint>

// Simple bit check function
inline constexpr bool check_bit (int input, int pos) {
  return (input & (1 << pos)) == (1 << pos);
}

// Simple test against floating zero
inline bool check_zero(double input, double eps = EPS) {
  return ( (input < eps) && (input > -eps) );
}


// Keeping this enum old-style, just for easier usage (not so
// verbose), such as boundary[LOWER_1] instead of
// boundary[(int)BoundaryPosition::LOWER_1]
// enum BoundaryPosition {
//   LOWER_1 = 0,
//   UPPER_1,
//   LOWER_2,
//   UPPER_2,
//   LOWER_3,
//   UPPER_3
// };         // ----- end of enum BoundaryPosition -----
using BoundaryPosition = int;
constexpr int LOWER_1 = 0;
constexpr int UPPER_1 = 1;
constexpr int LOWER_2 = 2;
constexpr int UPPER_2 = 3;
constexpr int LOWER_3 = 4;
constexpr int UPPER_3 = 5;

enum class FieldBCType {
  ROTATING_CONDUCTOR,
  DAMPING,
  COORDINATE,
};         // ----- end of enum FieldBCType -----

////////////////////////////////////////////////////////////////////////////////
///  Struct for storing the extent (3-dimensional size) of a multi-array.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// When computing curl, always Plus for B, and Minus for E
/// When computing divergence, Minus for B and Plus for E
/// When computing gradient, Plus for staggered scalar field and Minus
/// for unstaggered field
////////////////////////////////////////////////////////////////////////////////
enum class StaggerType : std::int8_t { STAGGER_MINUS, STAGGER_PLUS };

enum class FieldType : std::int8_t { ETYPE, BTYPE };

////////////////////////////////////////////////////////////////////////////////
///  Vectorized type as 3d vector
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct Vec3 {
  T x, y, z;

  typedef Vec3<T> self_type;
  typedef T elem_type;

  Vec3() : x(static_cast<T>(0)), y(static_cast<T>(0)), z(static_cast<T>(0)) {}
  Vec3(T xi, T yi, T zi) : x(xi), y(yi), z(zi) {}

  Vec3(const self_type& other) = default;
  Vec3(self_type&& other) = default;

  inline T& operator[](int idx) {
    if (idx == 0)
      return x;
    else if (idx == 1)
      return y;
    // else if (idx == 2)
    else
      return z;
    // else
    //   throw std::out_of_range("Index out of bound!");
  }

  inline const T& operator[](int idx) const {
    if (idx == 0)
      return x;
    else if (idx == 1)
      return y;
    // else if (idx == 2)
    else
      return z;
    // else
      // throw std::out_of_range("Index out of bound!");
  }

  inline self_type& operator= (const self_type& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  inline bool operator== (const self_type& other) const {
    return (x == other.x && y == other.y && z == other.z);
  }

  inline bool operator!= (const self_type& other) const {
    return (x != other.x || y != other.y || z != other.z);
  }

  // inline self_type& operator+=(const self_type& other) {
  //   x += other.x;
  //   y += other.y;
  //   z += other.z;
  //   return (*this);
  // }

  // inline self_type& operator-=(const self_type& other) {
  //   x -= other.x;
  //   y -= other.y;
  //   z -= other.z;
  //   return (*this);
  // }

  // inline self_type operator+(const self_type& other) const {
  //   self_type tmp = *this;
  //   tmp += other;
  //   return tmp;
  // }

  // inline self_type operator-(const self_type& other) const {
  //   self_type tmp = *this;
  //   tmp -= other;
  //   return tmp;
  // }


};

////////////////////////////////////////////////////////////////////////////////
///  Class to store a multi-dimensional size
////////////////////////////////////////////////////////////////////////////////
struct Extent : public Vec3<int>
{
  Extent() : Vec3(0, 1, 1) {}
  Extent(int w, int h = 1, int d = 1) : Vec3(w, h, d) {}
  Extent(const Vec3<int>& vec) : Vec3(vec) {}

  inline int& width() { return x; }
  inline const int& width() const { return x; }
  inline int& height() { return y; }
  inline const int& height() const { return y; }
  inline int& depth() { return z; }
  inline const int& depth() const { return z; }

  inline int size() const { return x * y * z; }

};

////////////////////////////////////////////////////////////////////////////////
///  Class to store a multi-dimensional index
////////////////////////////////////////////////////////////////////////////////
struct Index : public Vec3<int>
{
  Index() : Vec3(0, 0, 0) {}
  Index(int idx, const Extent& ext) {
    z = idx / (ext.width() * ext.height());
    int slice = idx % (ext.width() * ext.height());
    y = slice / ext.width();
    x = slice % ext.width();
  }
  Index(int xi, int yi, int zi) : Vec3(xi, yi, zi) {}
  Index(const Vec3<int>& vec) : Vec3(vec) {}

  inline int index(const Extent& ext) const {
    return x + y * ext.width() + z * ext.width() * ext.height();
  }

};

#endif  // ----- #ifndef _TYPES_H_  -----
