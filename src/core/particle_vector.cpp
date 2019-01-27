#include "core/particle_vector.hpp"
#include <stdexcept>
#include <algorithm>

namespace particle {
  template < typename T, std::size_t DPtc, typename state_t >
  vector<T, DPtc, state_t>::vector( std::size_t capacity )
    : _capacity(capacity), _size(0) {
    auto f =
      []( auto*& p, std::size_t n ) {
        p = new vec::remove_cvref_t<decltype(*p)> [n];
      };
    for ( auto*& ptr : _q ) f( ptr, capacity );
    for ( auto*& ptr : _p ) f( ptr, capacity );
    f( _state, capacity );
  }

  template < typename T, std::size_t DPtc, typename state_t >
  vector<T, DPtc, state_t>::~vector() {
    for ( auto*& ptr : _q ) delete [] ptr;
    for ( auto*& ptr : _p ) delete [] ptr;
    delete [] _state;
  }

  // real particles
  template < typename T, std::size_t DPtc, typename state_t >
  void vector<T, DPtc, state_t>::push_back( const Particle<T, DPtc, state_t>& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    vec::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
          q_arr[size] = q_ptc;
          p_arr[size] = p_ptc;
        }, _q, _p, ptc.q, ptc.p );
    _state[_size] = ptc.state;
    ++_size;
  }

  template < typename T, std::size_t DPtc, typename state_t >
  void vector<T, DPtc, state_t>::push_back( Particle<T, DPtc, state_t>&& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    vec::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, auto&& q_ptc, auto&& p_ptc ) {
          // simple types are better swapped than moved
          std::swap( q_arr[size], q_ptc );
          std::swap( p_arr[size], p_ptc );
        }, _q, _p, ptc.q, ptc.p );
    std::swap( _state[_size], ptc.state );
    ++_size;
  }

  // virtual particles
  template < typename T, std::size_t DPtc, typename state_t >
  void vector<T, DPtc, state_t>::push_back( const Particle< const T&, DPtc, const state_t& >& ptc ) {
    if ( _size == _capacity ) throw std::runtime_error( "size exceeds capacity" );
    vec::foreach<0,DPtc>
      ( [size=this->_size]( auto& q_arr, auto& p_arr, const auto& q_ptc, const auto& p_ptc ) {
          q_arr[size] = q_ptc;
          p_arr[size] = p_ptc;
        }, _q, _p, ptc.q, ptc.p );
    _state[_size] = ptc.state;
    ++_size;
  }


  template < typename T, std::size_t DPtc, typename state_t,
             typename PtcRef = vector<T, DPtc, state_t >::interator_t::reference >
  void swap( PtcRef a, PtcRef b ) noexcept {
    vec::swap( a.q, b.q );
    vec::swap( a.p, b.p );
    std::swap( a.state, b.state );
  }
}

namespace std {
  template < typename T, std::size_t DPtc, typename state_t,
             class Container = particle::vector<T,DPtc,state_t> >
  back_insert_iterator<Container>& back_insert_iterator<Container>::operator= ( const Particle< T, DPtc, state_t >& ptc ) {
    _c.push_back(ptc);
    return *this;
  }

  template < typename T, std::size_t DPtc, typename state_t,
             class Container = particle::vector<T,DPtc,state_t> >
  back_insert_iterator<Container>& back_insert_iterator<Container>::operator= ( const Particle< const T&, DPtc, const state_t& >& ptc ) {
    _c.push_back(ptc);
    return *this;
  }

  template < typename T, std::size_t DPtc, typename state_t,
             class Container = particle::vector<T,DPtc,state_t> >
  back_insert_iterator<Container>& back_insert_iterator<Container>::operator= ( Particle< T, DPtc, state_t >&& ptc ) {
    _c.push_back(std::move(ptc));
    return *this;
  }
}
