#include "testfw/testfw.hpp"
#include "particle/map.hpp"
#include "particle/array_impl.hpp"

using namespace particle;

template < typename T >
struct S {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = apt::copy_cvref_t<T, long long>;
};

template < typename T >
struct M {};

SCENARIO("Test structured binding syntax", "[particle]") {
  map<array<float,S>> ptcs;
  ptcs.insert( species::electron, array<float,S>() );
  ptcs.insert( species::ion, array<float,S>() );
  std::unordered_map<int,int> aaa;
  aaa.emplace( 56, 23 );
  aaa.emplace( 98, 42 );
  // for ( auto&[sp, ignore] : ptcs ) {
  //   std::cout << static_cast<int>(sp) << std::endl;
  //   ignore.resize(10);
  //   M<decltype(ignore)>::aa;
  // }
  for ( auto&[sp, ignore] : aaa ) {
    // std::cout << sp << std::endl;
    ignore = 656565;
    // M<decltype(ignore)>::aa;
  }
  std::cout << aaa[56] << ", " << aaa[98] << std::endl;
  // for ( auto&[sp, ignore] : aaa ) {
  //   std::cout << sp << std::endl;
  //   ignore = 10;
  //   M<decltype(ignore)>::aa;
  // }
}
