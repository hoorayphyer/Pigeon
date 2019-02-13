#include "particle/array.hpp"
#include <stdexcept>

namespace particle {
  template < typename T, int DPtc, typename state_t >
  array<T, DPtc, state_t >::array( std::size_t capacity )
    : _capacity(capacity), _size(0) {

    for ( T*& ptr : _q ) ptr = new T[capacity];
    for ( T*& ptr : _p ) ptr = new T[capacity];
    _state = new state_t[capacity];
  }

  template < typename T, int DPtc, typename state_t >
  array<T, DPtc, state_t>::~array() {
    for ( auto*& ptr : _q ) delete [] ptr;
    for ( auto*& ptr : _p ) delete [] ptr;
    delete [] _state;
  }

  // real particles
  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::push_back( const Particle<T, DPtc, state_t>& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    apt::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
          q_arr[size] = q_ptc;
          p_arr[size] = p_ptc;
        }, _q, _p, ptc.q, ptc.p );
    _state[_size] = ptc.state;
    ++_size;
  }

  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::push_back( Particle<T, DPtc, state_t>&& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    apt::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, auto&& q_ptc, auto&& p_ptc ) {
          // simple types are better swapped than moved
          std::swap( q_arr[size], q_ptc );
          std::swap( p_arr[size], p_ptc );
        }, _q, _p, ptc.q, ptc.p );
    std::swap( _state[_size], ptc.state );
    ++_size;
  }

  // // virtual particles
  // template < typename T, int DPtc >
  // void array<T, DPtc>::push_back( const Particle< const T&, DPtc >& ptc ) {
  //   if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
  //   apt::foreach<0,DPtc>
  //     ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
  //         q_arr[size] = q_ptc;
  //         p_arr[size] = p_ptc;
  //       }, _q, _p, ptc.q, ptc.p );
  //   _state[_size] = ptc.state;
  //   ++_size;
  // }

  template < typename T, int DPtc, typename state_t >
  void array<T, DPtc, state_t>::erase( int from, int to ) {
    auto min = []( auto a, auto b ) noexcept {
                 return ( a < b ) ? a : b; };
    if ( from > to ) {
      std::swap(from, to);
      ++from, ++to; // correct the -clusiveness
    }
    from = min( from, _size );
    to = min( to, _size );
    for ( int i = from; i != to; ++i )
      (*this)[i].set(flag::empty);
  }

}
