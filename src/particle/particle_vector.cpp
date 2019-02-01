#include "particle/particle_vector.hpp"
#include <stdexcept>
#include "apt/algorithm.hpp"

namespace particle {
  template < typename T, std::size_t DPtc >
  vector<T, DPtc >::vector( std::size_t capacity )
    : _capacity(capacity), _size(0) {
    auto f =
      []( auto*& p, std::size_t n ) {
        p = new apt::remove_cvref_t<decltype(*p)> [n];
      };
    for ( auto*& ptr : _q ) f( ptr, capacity );
    for ( auto*& ptr : _p ) f( ptr, capacity );
    f( _state, capacity );
  }

  template < typename T, std::size_t DPtc >
  vector<T, DPtc>::~vector() {
    for ( auto*& ptr : _q ) delete [] ptr;
    for ( auto*& ptr : _p ) delete [] ptr;
    delete [] _state;
  }

  // real particles
  template < typename T, std::size_t DPtc >
  void vector<T, DPtc>::push_back( const Particle<T, DPtc>& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    apt::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
          q_arr[size] = q_ptc;
          p_arr[size] = p_ptc;
        }, _q, _p, ptc.q, ptc.p );
    _state[_size] = ptc.state;
    ++_size;
  }

  template < typename T, std::size_t DPtc >
  void vector<T, DPtc>::push_back( Particle<T, DPtc>&& ptc ) {
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

  // virtual particles
  template < typename T, std::size_t DPtc >
  void vector<T, DPtc>::push_back( const Particle< const T&, DPtc >& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    apt::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
          q_arr[size] = q_ptc;
          p_arr[size] = p_ptc;
        }, _q, _p, ptc.q, ptc.p );
    _state[_size] = ptc.state;
    ++_size;
  }


  template < typename T, std::size_t DPtc,
             typename PtcRef = typename vector<T, DPtc >::interator_t::reference >
  void swap( PtcRef a, PtcRef b ) noexcept {
    apt::swap( a.q, b.q );
    apt::swap( a.p, b.p );
    std::swap( a.state, b.state );
  }
}
