#include "particle/array.hpp"

namespace particle {
  template < typename T, template < typename > class Specs >
  void array<T, Specs>::erase( unsigned int from, unsigned int to ) {
    auto min = []( auto a, auto b ) noexcept {
                 return ( a < b ) ? a : b; };
    if ( from > to ) {
      std::swap(from, to);
      ++from, ++to; // correct the -clusiveness
    }
    from = min( from, size() );
    to = min( to, size() );
    for ( int i = from; i != to; ++i )
      (*this)[i].reset(flag::exist);
  }

  template < typename T, template < typename > class Specs >
  void array<T, Specs>::resize(std::size_t size) {
    for ( int i = 0; i < DPtc; ++i ) {
      _q[i].resize( size, {} );
      _p[i].resize( size, {} );
    }
    _frac.resize(size, {});
    _state.resize(size, {});
  }
}
