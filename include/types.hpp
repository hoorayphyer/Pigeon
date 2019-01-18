#ifndef  _TYPES_HPP_
#define  _TYPES_HPP_

using Real = double; // type for a real number

constexpr Real
PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

struct ct_string { // compile time string
  char data[50];
};

constexpr bool strings_equal(char const * a, char const * b) {
  return *a == *b && ( *a == '\0' || strings_equal(a + 1, b + 1) );
}

constexpr bool operator== ( const ct_string& lhs, const ct_string& rhs ) {
  return strings_equal(lhs.data, rhs.data);
}

constexpr bool operator== ( const ct_string& lhs, const char* rhs ) {
  return strings_equal(lhs.data, rhs);
}

constexpr bool operator== ( const char* lhs, const ct_string& rhs ) {
  return strings_equal(lhs, rhs.data);
}

struct Species {
  ct_string str {"Unspecified"};
  int charge = 0; // in terms of unit charge
  int mass = 1; // in terms of unit mass
  bool is_radiative = false;

  // implicit conversion
  constexpr operator const char*() const {
    return str.data;
  }

};


constexpr Species Electron = { "Electron", -1, 1, true };
constexpr Species Positron = { "Positron", 1, 1, true };
constexpr Species Photon = { "Photon", 0, 0, false };

#include <unordered_map>
template < typename T >
using species_map = std::unordered_map<Species, T, std::hash<const char*>, std::equal_to<const char*>  >;


#endif
