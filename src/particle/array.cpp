#include "particle/array.hpp"
#include <stdexcept>

namespace particle {
  template < typename T, int DPtc, typename state_t >
  template < typename Ptc >
  void array<T, DPtc, state_t>::push_back( const PtcExpression<Ptc>& ptc ) {
    auto f = []( auto& arr, const auto& x ) { arr.push_back(x); };
    apt::foreach<0,DPtc> ( f, _q, ptc.q() );
    apt::foreach<0,DPtc> ( f, _p, ptc.p() );
    f( _state, ptc.state() );
  }

  template < typename T, int DPtc, typename state_t >
  template < typename Ptc >
  void array<T, DPtc, state_t>::push_back( PtcExpression<Ptc>&& ptc ) {
    auto f = []( auto& arr, auto&& x ) { arr.push_back(std::move(x)); };
    apt::foreach<0,DPtc> ( f, _q, ptc.q() );
    apt::foreach<0,DPtc> ( f, _p, ptc.p() );
    f( _state, ptc.state() );
  }

  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::erase( int from, int to ) {
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
    for ( auto& x : _q ) x.resize( size, {} );
    for ( auto& x : _p ) x.resize( size, {} );
    _state.resize(size, {});
    // flag padded particles as empty
    for ( auto i = old_size; i < size; ++i ) {
      (*this)[i].state().set(flag::empty);
    }

  }
}
