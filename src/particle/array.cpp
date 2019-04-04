#include "particle/array.hpp"

namespace particle {
  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::erase( unsigned int from, unsigned int to ) {
    auto min = []( auto a, auto b ) noexcept {
                 return ( a < b ) ? a : b; };
    if ( from > to ) {
      std::swap(from, to);
      ++from, ++to; // correct the -clusiveness
    }
    from = min( from, size() );
    to = min( to, size() );
    for ( int i = from; i != to; ++i )
      (*this)[i].set(flag::empty);
  }

  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::resize(std::size_t size) {
    const auto old_size = _state.size();
    for ( int i = 0; i < DPtc; ++i ) {
      _q[i].resize( size, {} );
      _p[i].resize( size, {} );
    }
    _state.resize(size, {});
    // flag padded particles as empty
    for ( auto i = old_size; i < size; ++i ) {
      (*this)[i].set(flag::empty);
    }

  }
}
